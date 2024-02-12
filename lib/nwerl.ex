defmodule NWerl do
  @moduledoc """
  Documentation for `NWerl`, a rudimentary Elixir module for sequence alignment using the
  [Needleman-Wunsch algorithm](https://en.wikipedia.org/w/index.php?title=Needleman%E2%80%93Wunsch_algorithm&useskin=vector).
  """

  @typedoc """
  A match scorer must take two terms and return an integer reflecting the match score for those terms.
  """
  @type match_scorer :: (term, term -> integer)

  @doc """
  Align two sequences using the
  [Needleman-Wunsch algorithm](https://en.wikipedia.org/w/index.php?title=Needleman%E2%80%93Wunsch_algorithm&useskin=vector)
  for sequence alignment.


  `NWerl.simple_match_scorer/3` will be used to score matches by default, but any
  function conforming to `t:NWerl.match_scorer/0` can be used.

  ## Examples

      iex> NWerl.align(["A", "B", "C"], ["B", "C", "D", "E"])
      [
        [
          {"A", nil},
          {"B", "B"},
          {"C", "C"},
          {nil, "D"},
          {nil, "E"}
        ]
      ]

  """
  @spec align([...], [...], match_scorer, integer, integer) :: [{any, any}, ...]
  def align(
        first_sequence,
        second_sequence,
        match_scorer \\ simple_match_scorer(),
        gap_open_penalty \\ -5,
        gap_extend_penalty \\ -1
      ) do
    scores = :ets.new(:scores, [:ordered_set])
    traceback = :ets.new(:traceback, [:bag])

    rows = Enum.count(first_sequence) + 1
    cols = Enum.count(second_sequence) + 1

    :ets.insert(scores, {{0, 0}, 0})
    :ets.insert(traceback, {{0, 0}, :DONE})

    1..rows
    |> Enum.each(fn i ->
      :ets.insert(scores, {{i, 0}, gap_open_penalty + gap_extend_penalty * i})
      :ets.insert(traceback, {{i, 0}, :UP})
    end)

    1..cols
    |> Enum.each(fn j ->
      :ets.insert(scores, {{0, j}, gap_open_penalty + gap_extend_penalty * j})
      :ets.insert(traceback, {{0, j}, :LEFT})
    end)

    # IO.puts(debug_print(scores, rows, cols, 0))

    1..rows
    |> Enum.each(fn i ->
      1..cols
      |> Enum.each(fn j ->
        potential_states =
          [:UP_LEFT, :UP, :LEFT]
          |> Enum.map(fn potential_state ->
            prior_score =
              case traverse_table(scores, {i, j}, potential_state) do
                [score] -> score
                _ -> 0
              end

            prior_state = traverse_table(traceback, {i, j}, potential_state)

            score =
              prior_score +
                cond do
                  potential_state == :UP_LEFT ->
                    match_scorer.(Enum.at(first_sequence, i - 1), Enum.at(second_sequence, j - 1))

                  (potential_state == :UP and prior_state == [:UP]) or
                      (potential_state == :LEFT and prior_state == [:LEFT]) ->
                    gap_extend_penalty

                  true ->
                    gap_open_penalty
                end

            {potential_state, score}
          end)

        {_, high_score} = potential_states |> Enum.max_by(fn {_state, score} -> score end)

        potential_states
        |> Enum.filter(fn {_state, score} -> score == high_score end)
        |> Enum.map(fn {state, _score} -> state end)
        |> Enum.each(fn state ->
          :ets.insert(traceback, {{i, j}, state})
        end)

        :ets.insert(scores, {{i, j}, high_score})
      end)
    end)

    # IO.puts(debug_print(scores, rows, cols, 0))
    # IO.puts(debug_print(traceback, rows, cols, ""))

    bottom_right = {rows - 1, cols - 1}

    possible_alignments = do_traceback(bottom_right, scores, traceback, [])

    # IO.inspect(possible_alignments)

    out =
      Enum.map(
        possible_alignments,
        fn possible_alignment ->
          possible_alignment
          |> Enum.map_reduce(bottom_right, fn e, acc ->
            {i, j} = acc

            {case e do
               :LEFT -> {nil, j - 1}
               :UP -> {i - 1, nil}
               :UP_LEFT -> {i - 1, j - 1}
             end, traverse(acc, e)}
          end)
          |> then(fn {list, _acc} -> list end)
          |> Enum.reverse()
          |> Enum.map(fn {a, b} ->
            {
              case a do
                nil -> nil
                i -> Enum.at(first_sequence, i)
              end,
              case b do
                nil -> nil
                i -> Enum.at(second_sequence, i)
              end
            }
          end)
        end
      )

    :ets.delete(scores)
    :ets.delete(traceback)

    out
  end

  @doc """
  Score a match between two elements using the provided match and mismatch scores.

  A comparator can optionally be provided; otherwise `Kernel.==/2` will be used to
  test for equality.
  """
  @spec simple_match_scorer(integer, integer, (term, term -> boolean)) :: match_scorer
  def simple_match_scorer(match_score \\ 10, mismatch_score \\ -100, comparator \\ &==/2) do
    fn first_thing, second_thing ->
      if comparator.(first_thing, second_thing) do
        match_score
      else
        mismatch_score
      end
    end
  end

  # defp debug_print(tab, rows, cols, default) do
  #   Enum.map(
  #     0..rows,
  #     fn i ->
  #       Enum.map(0..cols, fn j ->
  #         case :ets.lookup(tab, {i, j}) do
  #           [] -> default
  #           x when is_list(x) -> Enum.map(x, fn {_k, v} -> inspect(v) end) |> Enum.join(" ")
  #         end
  #       end)
  #       |> Enum.join("\t\t")
  #     end
  #   )
  #   |> Enum.join("\n")
  # end

  defp traverse({i, j}, :UP), do: {i - 1, j}
  defp traverse({i, j}, :LEFT), do: {i, j - 1}
  defp traverse({i, j}, :UP_LEFT), do: {i - 1, j - 1}

  defp traverse_table(table, starting_point, direction) do
    new_key = traverse(starting_point, direction)

    :ets.lookup(table, new_key)
    |> Enum.map(fn {_k, v} -> v end)
  end

  defp do_traceback(starting_point, scores, states, moves) do
    case :ets.lookup(states, starting_point) do
      [{_k, :DONE}] ->
        [moves |> Enum.reverse()]

      v ->
        v
        |> Enum.map(fn {_k, v} -> v end)
        |> Enum.flat_map(fn state ->
          next_point = traverse(starting_point, state)

          do_traceback(next_point, scores, states, [state | moves])
        end)
    end
  end
end
