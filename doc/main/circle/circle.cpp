#include "../../../src/cpp/solve/solve.cpp"

int main() {
  std::vector<std::string> matrixNames = {"_circleHandmade_3_2_6"};

  std::vector<Method> methods = {CI_FR, CN_FR, CI_L_BFGS, CN_L_BFGS};

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
      // https://stackoverflow.com/questions/8049556/what-s-the-difference-between-srand1-and-srand0
      for (int seed = 1; seed <= 10; seed++) {
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
        if (seed == 1)
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
