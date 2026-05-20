#pragma once

#include "../include/LBFGS.h"

#include <Eigen/Core>
#include <iostream>

#include "../util/problem.hpp"
#include "../util/timer.hpp"

template <typename ProblemT>
struct LBFGSFunction
{
  const ProblemT &problem;
  explicit LBFGSFunction(const ProblemT &problem) : problem(problem) {}

  double operator()(const Eigen::VectorXd &x, Eigen::VectorXd &grad) const
  {
    return static_cast<double>(problem.calc_score_and_grad(x, grad));
  }
};

template <typename ProblemT>
void solve_LBFGS(const ProblemT &problem, std::vector<Eigen::VectorXd> &positions,
                 std::vector<std::pair<double, double>> &hist, Timer &timer,
                 const int MAX_ITER)
{
  timer.start();
  LBFGSpp::LBFGSParam<double> param;
  param.m = 10;
  param.max_iterations = MAX_ITER;
  if (problem.modelType != ForceModelType::FR)
    param.max_iterations = 2 * MAX_ITER; // since non-FR models are harder to optimize, we allow more iterations
  param.epsilon = 1e-4;                  // to avoid line search failure
  param.epsilon_rel = 1e-4;              // to avoid line search failure
  LBFGSpp::LBFGSSolver<double> solver(param);

  LBFGSFunction<ProblemT> fun(problem);
  assert(!positions.empty());
  Eigen::VectorXd x = positions.back();
  double fx;

  try
  {
    solver.minimize(fun, x, fx, hist, positions, timer);
  }
  catch (const std::exception &e)
  {
    std::cerr << "\033[1;31mException: " << e.what() << "\033[0m" << std::endl;
  }

  timer.stop();
}
