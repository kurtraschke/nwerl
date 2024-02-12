defmodule NWerlTest do
  use ExUnit.Case
  doctest NWerl

  test "can align sequence" do
    assert NWerl.align(
      ["A", "B", "C", "D", "E", "F", "G"],
      ["B", "C", "D", "X", "G", "H"]
    ) ==
    [
      [
        {"A", nil},
        {"B", "B"},
        {"C", "C"},
        {"D", "D"},
        {"E", nil},
        {"F", nil},
        {nil, "X"},
        {"G", "G"},
        {nil, "H"}
      ]
    ]

  end
end
