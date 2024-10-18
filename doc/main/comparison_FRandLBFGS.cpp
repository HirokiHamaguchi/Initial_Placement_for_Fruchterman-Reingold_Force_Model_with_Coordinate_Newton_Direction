#include "../../src/cpp/solve/solve.cpp"

int main() {
  // define a problem (circle graph)
  size_t n = 6, m = 6;
  double k = 1;
  std::vector<size_t> row = {0, 1, 2, 3, 4, 0};
  std::vector<size_t> col = {1, 2, 3, 4, 5, 5};
  std::vector<double> data = {1, 1, 1, 1, 1, 1};
  Problem problem(n, m, k, row, col, data);

  // compute the LBFGS step from t=0
  std::srand(30);
  Eigen::VectorXf posInit = Eigen::VectorXf::Random(2 * problem.n);
  std::vector<Eigen::VectorXf> positions = {posInit};
  std::vector<std::pair<double, double>> _hist;
  Timer _timer;
  solve_LBFGS<FunctionFR>(problem, positions, _hist, _timer, 10);

  int t = 6;

  // output position at t
  std::cout << "posInit=[";
  for (size_t i = 0; i < problem.n; ++i) {
    std::cout << "(" << positions[t][2 * i] << "," << positions[t][2 * i + 1] << "),";
  }
  std::cout << "]" << std::endl;

  // output position at t+1
  std::cout << "posLBFGS=[";
  for (size_t i = 0; i < problem.n; ++i) {
    std::cout << "(" << positions[t + 1][2 * i] << "," << positions[t + 1][2 * i + 1]
              << "),";
  }
  std::cout << "]" << std::endl;

  // compute the FR step from t
  std::vector<Eigen::VectorXf> positions2 = {positions[t]};
  std::vector<std::pair<double, double>> _hist2;
  solve_FR(problem, positions2, _hist2, _timer, 10);
  assert(positions2[0] == positions[t]);

  // output position at t+1 (i.e. positions2[1])
  std::cout << "posFR=[";
  for (size_t i = 0; i < problem.n; ++i) {
    std::cout << "(" << positions2[1][2 * i] << "," << positions2[1][2 * i + 1] << "),";
  }
  std::cout << "]" << std::endl;

  // L-BFGS output is better than FR output
  dbg(problem.calcScore(positions[t]));
  dbg(problem.calcScore(positions[t + 1]));
  dbg(problem.calcScore(positions2[0]));
  dbg(problem.calcScore(positions2[1]));

  return 0;
}
