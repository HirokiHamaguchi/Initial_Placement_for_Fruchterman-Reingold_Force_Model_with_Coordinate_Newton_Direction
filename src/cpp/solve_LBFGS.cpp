#pragma once

#include <LBFGS.h>

#include <Eigen/Core>
#include <iostream>

#include "util/function.hpp"
#include "util/problem.hpp"

Eigen::VectorXf solve_LBFGS(const Problem& problem, const Eigen::VectorXf& x_init) {
  LBFGSpp::LBFGSParam<float> param;
  param.m = 10;
  param.max_iterations = 100;
  LBFGSpp::LBFGSSolver<float> solver(param);

  Function fun(problem);
  Eigen::VectorXf x = x_init;
  float fx;

  auto t0 = std::chrono::high_resolution_clock::now();
  int niter = solver.minimize(fun, x, fx);
  auto t1 = std::chrono::high_resolution_clock::now();

  std::cout << "Elapsed time: "
            << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
            << " seconds\n"
            << " iterations : " << niter << "\n"
            << "f(x) = " << fx << "\n"
            << "||grad|| = " << solver.final_grad_norm() << std::endl;

  return x;
}
