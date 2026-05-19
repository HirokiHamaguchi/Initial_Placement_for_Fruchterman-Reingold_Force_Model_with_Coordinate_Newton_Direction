#pragma once

#include <Eigen/Core>
#include <vector>
#include <cmath>

/**
 * Abstract base class for graph drawing force models (FR, HC, Eades)
 * Defines the interface for energy computation, gradients, Hessians, and scaling
 */
class ForceModel
{
public:
  virtual ~ForceModel() = default;

  // Network structure accessors
  size_t getN() const { return n; }
  size_t getM() const { return m; }

  /**
   * Calculate total energy (attractive + repulsive)
   * Used for scoring during optimization
   * @param position Current node positions (2*n vector)
   * @param includeRepulsive If false, compute only attractive energy
   * @return Energy value
   */
  virtual double calcScore(const Eigen::VectorXf &position,
                           bool includeRepulsive = true) const = 0;

  /**
   * Calculate energy and gradient simultaneously (efficient for LBFGS)
   * @param x Current node positions (2*n vector)
   * @param grad Output: gradient vector (2*n)
   * @return Energy value
   */
  virtual double calc_score_and_grad(const Eigen::VectorXf &x,
                                     Eigen::VectorXf &grad) const = 0;

  /**
   * Calculate gradient and Hessian for coordinate update
   * Computes local Hessian for each coordinate and applies it in place
   * @param dist Distance between two nodes
   * @param dx, dy Displacement vector components
   * @param w Edge weight (or target distance)
   * @param gx, gy Output: gradient components
   * @param hxx, hxy, hyy Output: Hessian elements (upper triangle)
   */
  virtual void calc_grad_hess(float dist, float dx, float dy, float w,
                              float &gx, float &gy, float &hxx, float &hxy,
                              float &hyy) const = 0;

  /**
   * Optimal global scaling of positions
   * Minimizes energy over scale factor s: phi(s) = Phi(s*X)
   * @param position Node positions (modified in place)
   */
  virtual void optimalScaling(Eigen::VectorXf &position) const = 0;

protected:
  size_t n;                                                // number of vertices
  size_t m;                                                // number of edges
  std::vector<size_t> row;                                 // edge u
  std::vector<size_t> col;                                 // edge v
  std::vector<double> data;                                // edge weight (= target distance in HC/Eades)
  std::vector<std::vector<std::pair<size_t, double>>> adj; // adjacency list

  // Helper: build adjacency list from row/col/data
  void makeAdj()
  {
    adj.resize(n);
    for (size_t i = 0; i < m; ++i)
    {
      adj[row[i]].emplace_back(col[i], data[i]);
      adj[col[i]].emplace_back(row[i], data[i]);
    }
  }
};
