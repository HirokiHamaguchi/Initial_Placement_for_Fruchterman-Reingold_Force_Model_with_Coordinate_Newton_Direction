#include "solve_init.cpp"
#include "util/problem.hpp"

int main() {
  std::vector<size_t> row = {0, 0};
  std::vector<size_t> col = {1, 2};
  std::vector<double> data = {1.0, 1.0};
  Problem problem(3, 2, 1.0, row, col, data);

  solve_init(problem);

  return 0;
}