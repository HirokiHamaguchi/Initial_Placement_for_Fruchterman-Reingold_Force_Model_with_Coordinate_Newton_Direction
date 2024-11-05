#include "../../../src/cpp/solve/solve.cpp"

int main() {
  std::string matrixNamesPath = "matrixNames.txt";
  std::ifstream matrixNamesFile(matrixNamesPath);
  assert(matrixNamesFile.is_open());
  std::vector<std::string> matrixNames;
  std::string matrixName;
  while (std::getline(matrixNamesFile, matrixName)) matrixNames.push_back(matrixName);
  matrixNamesFile.close();

  std::vector<Method> methods = {CI_FR, CN_FR, CI_L_BFGS, CN_L_BFGS};

  for (auto& default_ITER : {15, 50}) {
    std::string histStr = std::to_string(matrixNames.size()) + " " +
                          std::to_string(methods.size()) + "\n";

    for (std::string matrixName : matrixNames) {
      Problem problem(matrixName);
      histStr += matrixName + "\n";
      std::cout << matrixName << std::endl;
      for (auto& method : methods) {
        int MAX_ITER =
            (method == CN_FR || method == CN_L_BFGS) ? default_ITER - 3 : default_ITER;

        histStr += MethodStr[method] + "\n";
        std::cout << MethodStr[method] << std::endl;
        // Solve by each method
        // https://stackoverflow.com/questions/8049556/what-s-the-difference-between-srand1-and-srand0
        for (int seed = 1; seed <= 10; seed++) {
          auto [hist, positions] = solve(method, problem, false, seed, MAX_ITER);
          // Output
          assert(!hist.empty());
          double score = hist.back().first;
          histStr += "Score: " + std::to_string(score) + "\n";
          std::cout << "Score: " << score << std::endl;
          if (seed == 1)
            problem.printOutput(positions, "out/" + problem.matrixName + "_" +
                                               MethodStr[method] + "_" +
                                               std::to_string(default_ITER) + ".out");
        }
      }
    }

    std::string fileName = "hist_" + std::to_string(default_ITER) + ".txt";
    auto [histPath, fileForHist] = openFile(fileName);
    fileForHist << histStr;
    fileForHist.close();
    std::cout << "Hist path: " << histPath << std::endl;
  }

  return 0;
}
