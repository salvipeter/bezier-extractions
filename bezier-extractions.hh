#pragma once

#include <vector>

#include <Eigen/Core>

std::vector<Eigen::MatrixXd> bezierExtractionMatrices(size_t p, const std::vector<double> &knots);
