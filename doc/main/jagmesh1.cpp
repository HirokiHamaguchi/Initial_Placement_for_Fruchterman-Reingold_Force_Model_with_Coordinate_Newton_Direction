#include "../../src/cpp/solve/solve.cpp"

int main() {
  // define a problem
  Problem problem("jagmesh1");

  // set seed and make random initial positions
  int seed = 0;
  std::srand(seed);
  Eigen::VectorXf posInit = Eigen::VectorXf::Random(2 * problem.n);
  std::vector<Eigen::VectorXf> positions;
  std::vector<std::pair<double, double>> _hist;
  Timer _timer;

  std::vector<Eigen::VectorXf> finalOutputs;

  const int MAX_ITER = 50;

  // FR
  positions = {posInit};
  _hist.clear();
  solve_FR(problem, positions, _hist, _timer, MAX_ITER);
  dbg(problem.calcScore(positions.back()));
  finalOutputs.push_back(positions.back());

  // L_BFGS
  positions = {posInit};
  _hist.clear();
  solve_LBFGS<FunctionFR>(problem, positions, _hist, _timer, MAX_ITER);
  dbg(problem.calcScore(positions.back()));
  finalOutputs.push_back(positions.back());

  // RS
  positions = {};
  _hist.clear();
  solve_init(problem, false, seed, positions, _hist, _timer);
  dbg(problem.calcScore(positions.back()));
  finalOutputs.push_back(positions.back());

  // RS_L_BFGS
  positions = {positions.back()};
  _hist.clear();
  solve_LBFGS<FunctionFR>(problem, positions, _hist, _timer, MAX_ITER);
  dbg(problem.calcScore(positions.back()));
  finalOutputs.push_back(positions.back());

  problem.printOutput(finalOutputs, "doc/main/jagmesh1.txt");

  return 0;
}
