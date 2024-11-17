#include "../../../src/cpp/solve/solve.cpp"

int main() {
  std::vector<Method> methods = {SA_L_BFGS, CN_L_BFGS};
  int default_ITER = 10;

  size_t n = 100, m = 1000;
  double k = 1.0 / std::sqrt(n);
  std::vector<size_t> row, col;
  std::vector<double> data;
  std::mt19937 mt(1);
  std::uniform_int_distribution<int> dist(0, n - 1);
  while (true) {
    int i = dist(mt), j = dist(mt);
    if (i == j) continue;
    if (i > j) std::swap(i, j);
    row.push_back(i);
    col.push_back(j);
    data.push_back(i / 33 == j / 33 ? 1.0 : 0.1);
    if (row.size() == m) break;
  }

  Problem problem(n, m, k, row, col, data);
  problem.matrixName = "dense";

  std::string histStr = "1 " + std::to_string(methods.size()) + "\n";
  histStr += "dense\n";

  for (auto& method : methods) {
    int MAX_ITER = default_ITER;
    histStr += MethodStr[method] + "\n";
    std::cout << MethodStr[method] << std::endl;
    int seed = 1;
    auto [hist, positions] = solve(method, problem, false, seed, MAX_ITER);
    assert(!hist.empty());
    double score = hist.back().first;
    histStr += "Score: " + std::to_string(score) + "\n";
    std::cout << "Score: " << score << std::endl;
    problem.printOutput(positions, "out/" + problem.matrixName + "_" +
                                       MethodStr[method] + "_" +
                                       std::to_string(default_ITER) + ".out");
  }

  std::string fileName = "hist2_" + std::to_string(default_ITER) + ".txt";
  auto [histPath, fileForHist] = openFile(fileName);
  fileForHist << histStr;
  fileForHist.close();
  std::cout << "Hist path: " << histPath << std::endl;

  return 0;
}
