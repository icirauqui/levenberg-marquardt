#include <iostream>
#include <random>
#include "mylib/mylib.hpp"

#include <cmath>
#include "eigen3/Eigen/Dense"


void LevenbergMarquardt(const Eigen::MatrixXd& data,
                        const Eigen::VectorXd& initial_params,
                        const int max_iterations,
                        const double epsilon,
                        double lambda,
                        Eigen::VectorXd& final_params) {

  std::cout << "LevenbergMarquardt " << std::endl;

  const int num_params = initial_params.size();
  Eigen::VectorXd params = initial_params;
  double error = std::numeric_limits<double>::max();
  int iteration = 0;

  while (error > epsilon && iteration < max_iterations) {
    Eigen::MatrixXd J(data.rows(), num_params);
    Eigen::VectorXd r(data.rows());

    for (int i = 0; i < data.rows(); i++) {
      const double x = data(i,0);
      const double y = data(i,1);
      const double f = std::exp(-params(0) * x) * std::cos(params(1) * x);
      
      r(i) = y - f;

      J(i, 0) = x * f;
      J(i, 1) = -x * std::exp(-params(0) * x) * std::sin(params(1) * x);
    }

    Eigen::MatrixXd A = J.transpose() * J + lambda * Eigen::MatrixXd::Identity(num_params, num_params);
    Eigen::VectorXd b = J.transpose() * r;

    Eigen::VectorXd delta = A.ldlt().solve(b);
    Eigen::VectorXd new_params = params + delta;

    const double new_error = (r.array() * r.array()).sum();

    if (new_error < error) {
      lambda /= 10.0;
      error = new_error;
      params = new_params;
    } else {
      lambda *= 10.0;
    }

    iteration++;
  }

  final_params = params;
}


int main() {
  std::cout << "Hello, World!" << std::endl;

  const int num_data_points = 100;
  const int max_iterations = 1000;
  const double lambda = 0.01;
  const double epsilon = 0.000001;
  Eigen::VectorXd initial_params(2);
  initial_params << 0.5, 1.0;
  Eigen::VectorXd final_params;

  // Generate data
  std::default_random_engine gen;
  std::normal_distribution<double> dist(0.0, 0.1);
  Eigen::MatrixXd data(num_data_points, 2);

  for (int i = 0; i < num_data_points; ++i) {
    data(i,0) = i / static_cast<double>(num_data_points - 1);
    data(i,1) = std::exp(-initial_params(0) * data(i,0)) * std::cos(initial_params(1) * data(i,0)) + dist(gen);
  }

  // Run Levenberg-Marquardt
  LevenbergMarquardt(data, initial_params, max_iterations, epsilon, lambda, final_params);

  // Print results
  std::cout << "Initial parameters: " << initial_params.transpose() << std::endl;
  std::cout << "Final parameters: " << final_params.transpose() << std::endl;
  
  //std::cout << "Sum of squared residuals: " << (data.col(1) - std::exp(-final_params(0) * data.col(0).array()) * std::cos(final_params(1) * data.col(0).array())).squaredNorm() << std::endl;
    


  int a = mylib();
  return 0;
}
