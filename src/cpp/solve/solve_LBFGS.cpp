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
  param.epsilon_rel = 1e-3;
  LBFGSpp::LBFGSSolver<float> solver(param);

  MyFunction fun(problem);
  assert(!positions.empty());
  Eigen::VectorXf x = positions.back();
  float fx;

  int niter = solver.minimize(fun, x, fx, hist, positions, timer);
  timer.stop();

  dbg(niter, fx, solver.final_grad_norm());
}
