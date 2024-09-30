#pragma once

#include <LBFGS.h>

#include <Eigen/Core>
#include <iostream>

#include "util/function.hpp"
#include "util/problem.hpp"

std::vector<Eigen::VectorXf> solve_LBFGS(const Problem& problem,
                                         const Eigen::VectorXf& x_init,
                                         const bool measureTime) {
  LBFGSpp::LBFGSParam<float> param;
  param.m = 10;
  param.max_iterations = 100;
  param.epsilon_rel = 1e-3;
  LBFGSpp::LBFGSSolver<float> solver(param);

  Function fun(problem);
  Eigen::VectorXf x = x_init;
  float fx;

  auto [niter, positions] = solver.minimize(fun, x, fx, measureTime);
  assert(!measureTime || int(positions.size()) <= 1);
  dbg(niter, fx, solver.final_grad_norm());

  return positions;
}
