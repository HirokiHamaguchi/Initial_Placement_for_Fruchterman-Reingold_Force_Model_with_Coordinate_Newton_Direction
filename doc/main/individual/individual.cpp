#include "../../../src/cpp/solve/solve.cpp"

int main() {
  std::vector<std::string> matrixNames = {
      "cycle300", "jagmesh1",  // my selection
      "btree9",                // from networkx
      "1138_bus", "dwt_1005", "dwt_2680",
      "3elt",  // SGD
  };

  std::vector<Method> methods = {FR, RS_FR, L_BFGS, RS_L_BFGS};

  std::string histStr =
      std::to_string(matrixNames.size()) + " " + std::to_string(methods.size()) + "\n";

  const int MAX_ITER = 200;

  for (std::string matrixName : matrixNames) {
    Problem problem(matrixName);
    histStr += matrixName + "\n";
    std::cout << matrixName << std::endl;
    for (auto& method : methods) {
      histStr += MethodStr[method] + "\n";
      std::cout << MethodStr[method] << std::endl;
      // Solve by each method
      for (int seed = 0; seed < 3; seed++) {
        auto [hist, positions] = solve(method, problem, true, seed, MAX_ITER);
        // Output
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
        if (seed == 0)
          problem.printOutput(positions, "out/" + problem.matrixName + "_" +
                                             MethodStr[method] + ".out");
      }
    }
  }

  std::string fileName = "hist.txt";
  auto [histPath, fileForHist] = openFile(fileName);
  fileForHist << histStr;
  fileForHist.close();
  std::cout << "Hist path: " << histPath << std::endl;

  return 0;
}
