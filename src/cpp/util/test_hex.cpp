#include <cassert>
#include <iostream>
#include <vector>

#include "hex.hpp"

int main() {
  std::vector<Hex> a = {Hex(0, 0, 0),  Hex(0, -1, 1), Hex(0, -2, 2),
                        Hex(1, -3, 2), Hex(1, -4, 3), Hex(1, -5, 4)};
  std::vector<Hex> b = Hex::linedraw(Hex(0, 0, 0), Hex(1, -5, 4));

  for (auto& h : b) std::cerr << h << std::endl;

  assert(Hex::equal_hex_array(a, b));

  return 0;
}