#include "solve.cpp"

int main(int argc, char* argv[]) {
  // Read input
  Method method;
  Problem problem;
  bool measureTime;

  if (argc == 1) {
    method = RS_L_BFGS;
    problem = Problem("jagmesh1");
    measureTime = true;
  } else if (argc == 3) {
    std::string problemStr = argv[1];
    std::string measureTimeStr = argv[2];
    problem = Problem(problemStr);
    if (measureTimeStr != "true" && measureTimeStr != "false") {
      std::cerr << "Error: invalid measureTime\n";
      exit(1);
    }
    measureTime = (measureTimeStr == "true");
  } else {
    std::cerr << "Error: invalid number of arguments\n";
    exit(1);
  }

  if (!measureTime) {
    auto [hist, positions] = solve(method, problem, measureTime, 0);
    assert(!hist.empty());
    std::cout << "Elapsed time: " << hist.back().second << " seconds" << std::endl;
    std::cout << "Score: " << hist.back().first << std::endl;
    problem.printOutput(positions);
    return 0;
  }

  std::vector<std::string> matrixNames = {
      // "jagmesh1", "jagmesh2", "jagmesh3", "dwt_221", "dwt_1005", "dwt_1005",
      // "arc130",
      // "ash85",    "ash292",   "bcspwr08", "bp_800",  "can_715",  "ash85"
      "dwt_1005",
  };

  std::vector<Method> methods = {FR, RS_FR, L_BFGS, RS_L_BFGS};

  std::string histStr =
      std::to_string(matrixNames.size()) + " " + std::to_string(methods.size()) + "\n";

  for (std::string matrixName : matrixNames) {
    problem = Problem(matrixName);
    histStr += matrixName + "\n";
    std::cout << matrixName << std::endl;
    for (auto& method : methods) {
      histStr += MethodStr[method] + "\n";
      std::cout << MethodStr[method] << std::endl;
      // Solve by each method
      const int seed = 12345;
      auto [hist, positions] = solve(method, problem, measureTime, seed);
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
      problem.printOutput(positions);
    }
  }
  std::string fileName = "out/_hist.txt";
  auto [histPath, fileForHist] = openFile(fileName);
  fileForHist << histStr;
  fileForHist.close();
  std::cout << "Hist path: " << histPath << std::endl;

  return 0;
}
