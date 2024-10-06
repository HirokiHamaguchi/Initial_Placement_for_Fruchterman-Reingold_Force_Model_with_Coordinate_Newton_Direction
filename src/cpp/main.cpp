#include "include/LBFGS.h"
#include "solve_FR.cpp"
#include "solve_LBFGS.cpp"
#include "solve_init.cpp"
#include "solve_init2.cpp"
#include "util/function.hpp"
#include "util/problem.hpp"

enum Method {
  FR,         // Fruchterman-Reingold method
  L_BFGS,     // Limited-memory Broyden-Fletcher-Goldfarb-Shanno method
  RS_FR,      // Init:RS / Optimize:FR
  RS_L_BFGS,  // Init:RS / Optimize:L_BFGS
};

std::vector<Eigen::VectorXf> main_sub(const Method method, const Problem& problem,
                                      const bool measureTime, const int seed) {
  std::vector<Eigen::VectorXf> positions;
  if (method == RS_FR || method == RS_L_BFGS) {
    positions = solve_init(problem, measureTime, seed);
  } else {
    std::srand(0);
    Eigen::VectorXf position = Eigen::VectorXf::Random(2 * problem.n);
    for (int i = 0; i < position.size(); ++i) position[i] = std::abs(position[i]);
    positions.push_back(position);
  }
  std::vector<Eigen::VectorXf> positions2;
  if (method == L_BFGS || method == RS_L_BFGS) {
    positions2 = solve_LBFGS<FunctionFR>(problem, positions.back(), measureTime);
  } else if (method == FR || method == RS_FR) {
    positions2 = solve_FR(problem, positions.back(), measureTime);
  }
  positions.insert(positions.end(), positions2.begin(), positions2.end());
  return positions;
}

int main(int argc, char* argv[]) {
  // Read input
  Method method = RS_L_BFGS;  // only for !measureTime
  Problem problem;
  bool measureTime;
  if (argc == 1) {
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

  if (measureTime) {
    for (auto& method : {L_BFGS, RS_L_BFGS}) {
      // Solve by each method
      double time = 0.0;
      double score = 0.0;
      std::vector<Eigen::VectorXf> positions;
      const int nTrials = 1;
      for (int seed = 0; seed < nTrials; seed++) {
        auto t0 = std::chrono::high_resolution_clock::now();
        auto newPositions = main_sub(method, problem, measureTime, seed);
        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsedTime = std::chrono::duration<double>(t1 - t0).count();
        std::cout << "Elapsed time: " << elapsedTime << " seconds" << std::endl;
        time += elapsedTime, score += problem.calcScore(newPositions.back());
        if (positions.empty()) positions = newPositions;
      }
      time /= nTrials, score /= nTrials;
      // Output
      std::cout << "Average Elapsed time: " << time << " seconds" << std::endl;
      std::cout << "Average Score: " << score << std::endl;
      problem.printOutput(positions);
    }
  } else {
    // Solve by the specified method
    auto t0 = std::chrono::high_resolution_clock::now();
    std::vector<Eigen::VectorXf> positions = main_sub(method, problem, measureTime, 0);
    auto t1 = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(t1 - t0).count();
    // Output
    std::cout << "Average Elapsed time: " << time << " seconds" << std::endl;
    std::cout << "Final Score: " << problem.calcScore(positions.back()) << std::endl;
    problem.printOutput(positions);
  }

  return 0;
}
