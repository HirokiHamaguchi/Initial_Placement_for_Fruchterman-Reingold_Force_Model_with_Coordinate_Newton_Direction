#pragma once

#include <LBFGS.h>

#include <Eigen/Core>
#include <iostream>

#include "util/problem.hpp"

template <typename MyFunction>
std::vector<Eigen::VectorXf> solve_LBFGS(const Problem& problem,
                                         const Eigen::VectorXf& x_init,
                                         const bool measureTime) {
  LBFGSpp::LBFGSParam<float> param;
  param.m = 10;
  param.max_iterations = 50;
  param.epsilon_rel = 1e-3;
  LBFGSpp::LBFGSSolver<float> solver(param);

  MyFunction fun(problem);
  Eigen::VectorXf x = x_init;
  float fx;

  auto [niter, positions] = solver.minimize(fun, x, fx, measureTime);
  assert(!measureTime || int(positions.size()) <= 1);
  dbg(niter, fx, solver.final_grad_norm());

  return positions;
}
