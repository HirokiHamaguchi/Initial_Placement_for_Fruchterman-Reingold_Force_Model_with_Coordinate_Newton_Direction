#include "../../../src/cpp/solve/solve.cpp"

int main() {
  std::vector<std::string> matrixNames = {
      "cycle300", "jagmesh1", "btree9", "1138_bus", "dwt_1005", "dwt_2680", "3elt",
  };

  std::vector<Method> methods = {FR, CN_FR, L_BFGS, CN_L_BFGS};

  std::string histStr =
      std::to_string(matrixNames.size()) + " " + std::to_string(methods.size()) + "\n";

  const int MAX_ITER = 200;

  for (auto matrixName : matrixNames) {
    Problem problem(matrixName);
    histStr += matrixName + "\n";
    std::cout << matrixName << std::endl;
    for (auto& method : methods) {
      histStr += MethodStr[method] + "\n";
      std::cout << MethodStr[method] << std::endl;

      // Solve by each method
      // https://stackoverflow.com/questions/8049556/what-s-the-difference-between-srand1-and-srand0
      for (int seed = 1; seed <= 10; seed++) {
        auto [hist, positions] = solve(method, problem, true, seed, MAX_ITER);
        histStr += std::to_string(hist.size()) + "\n";
        for (auto [score, _] : hist) histStr += std::to_string(score) + " ";
        histStr += "\n";
        for (auto [_, time] : hist) histStr += std::to_string(time) + " ";
        histStr += "\n";
        auto [score, time] = hist.back();
        histStr += "Elapsed_time: " + std::to_string(time) + "\n";
        histStr += "Score: " + std::to_string(score) + "\n";
        std::cout << "Elapsed time: " << time << " seconds" << std::endl;
        std::cout << "Score: " << score << std::endl;
        // Output the result
        if (seed == 1)
          problem.printOutput(positions, "out/" + problem.matrixName + "_" +
                                             MethodStr[method] + ".out");
      }

      if (method != CN_L_BFGS) continue;

      // Compute the optimal solution with the seed 1
      const int MAX_ITER_2 = 500;
      auto [_, medPos] = solve(method, problem, true, 1, MAX_ITER_2);
      problem.printOutput(medPos, "out/opt_" + problem.matrixName + ".out");
    }
  }

  std::string fileName = "hist.txt";
  auto [histPath, fileForHist] = openFile(fileName);
  fileForHist << histStr;
  fileForHist.close();
  std::cout << "Hist path: " << histPath << std::endl;

  return 0;
}
