#pragma once

#include <LBFGS.h>

#include <Eigen/Core>
#include <iostream>

#include "../util/problem.hpp"
#include "../util/timer.hpp"

template <typename MyFunction>
void solve_LBFGS(const Problem& problem, std::vector<Eigen::VectorXf>& positions,
                 std::vector<std::pair<double, double>>& hist, Timer& timer,
                 const int MAX_ITER) {
  timer.start();
  LBFGSpp::LBFGSParam<float> param;
  param.m = 10;
  param.max_iterations = MAX_ITER;
  param.epsilon = 1e-4;      // to avoid line search failure
  param.epsilon_rel = 1e-4;  // to avoid line search failure
  LBFGSpp::LBFGSSolver<float> solver(param);

  MyFunction fun(problem);
  assert(!positions.empty());
  Eigen::VectorXf x = positions.back();
  float fx;

  try {
    solver.minimize(fun, x, fx, hist, positions, timer);
  } catch (const std::exception& e) {
    std::cerr << "\033[1;31mException: " << e.what() << "\033[0m" << std::endl;
  }

  timer.stop();
}
