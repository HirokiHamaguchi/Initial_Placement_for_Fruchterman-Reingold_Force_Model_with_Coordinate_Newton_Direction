#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

#include "../util/grid.hpp"
#include "../util/hex.hpp"
#include "../util/problem.hpp"
#include "../util/timer.hpp"

void addVis(const Grid &grid, std::vector<Eigen::VectorXd> &positions, int it,
            bool measureTime, const int ITERATIONS)
{
  if (measureTime)
    return;
  if (it % (ITERATIONS / 10) != 0)
    return;
  positions.push_back(grid.toPosition());
}

void solve_init(const Problem &problem, const bool measureTime, const int seed,
                std::vector<Eigen::VectorXd> &positions)
{
  Grid grid(problem.n, seed);
  auto firstPos = grid.toPosition();
  if (problem.modelType == ForceModelType::HC || problem.modelType == ForceModelType::Eades)
  {
    double initScaling = problem.optimalScaling(firstPos);
    grid.scale(initScaling);
  }

  // initialize random number generator
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> distVertex(0, problem.n - 1);
  std::uniform_real_distribution<double> distHexR(0, 1.0);
  std::uniform_real_distribution<double> distHexTheta(0, 2 * M_PI);

  // Randomness parameters
  const int ITERATIONS = 2 * problem.n * (problem.n * problem.n / problem.m);
  const double T0 = +1.0;
  const double T1 = +0.5;

  for (int it = 0; it < ITERATIONS; it++)
  {
    // randomly select a vertex
    const size_t i = distVertex(gen);
    const Hex hexI = grid.points[i];

    // calculate gradient and Hessian
    double gx = 0.0, gy = 0.0;
    double hxx = 0.0, hxy = 0.0, hyy = 0.0;
    for (auto [j, w] : problem.adj[i])
    {
      int dq = hexI.q - grid.points[j].q;
      int dr = hexI.r - grid.points[j].r;
      grid.calc_grad_hess(problem, dq, dr, w, gx, gy, hxx, hxy, hyy);
    }

    // if local minimum, continue
    if ((gx * gx + gy * gy) < 1e-9)
    {
      addVis(grid, positions, it, measureTime, ITERATIONS);
      continue;
    }

    // compute Newton's direction (Hess^{-1} @ (-grad))
    auto [dx, dy] = grid.computeDxDy(gx, gy, hxx, hxy, hyy);
    double norm_d = std::hypot(dx, dy);
    if (problem.modelType == ForceModelType::HC || problem.modelType == ForceModelType::Eades)
    {
      // cap the step size to avoid instability
      if (norm_d > 3.0 * grid.grid_size)
      {
        double scale = (3.0 * grid.grid_size) / norm_d;
        dx *= scale;
        dy *= scale;
      }
    }

    // select random neighbor of (x + dx, y + dy) with randomness
    const auto [x, y] = grid.hex2xy(hexI.q, hexI.r);
    double T = std::max((T0 + (T1 - T0) * it / ITERATIONS) * norm_d, 1.0 * grid.grid_size);
    double r = T * distHexR(gen), theta = distHexTheta(gen);
    double dxr = r * std::cos(theta), dyr = r * std::sin(theta);
    const Hex hexJ = grid.xy2hex(x + dx + dxr, y + dy + dyr);

    // swap position
    if (grid.isInside(hexJ) && hexI != hexJ)
      grid.swap(i, hexI, hexJ);

    addVis(grid, positions, it, measureTime, ITERATIONS);

    // assert(grid.isCorrectState()); // ! For debug
  }

  auto finalPos = grid.toPosition();
  problem.optimalScaling(finalPos);
  positions.push_back(finalPos);
}
