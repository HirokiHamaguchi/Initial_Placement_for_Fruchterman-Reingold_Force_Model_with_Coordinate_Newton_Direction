#include "solve/solve.cpp"

int main() {
  Method method = FR;
  Problem problem("jagmesh1");
  bool measureTime = false;

  auto [hist, positions] = solve(method, problem, measureTime, 0, 100);
  assert(!hist.empty());
  std::cout << "Elapsed time: " << hist.back().second << " seconds" << std::endl;
  std::cout << "Score: " << hist.back().first << std::endl;
  problem.printOutput(positions,
                      "out/" + problem.matrixName + "_" + MethodStr[method] + ".out");

  return 0;
}
