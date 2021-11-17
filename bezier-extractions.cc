#include "bezier-extractions.hh"

using Eigen::MatrixXd;

// As described in
//   M.J. Borden, M.A. Scott, J.A. Evans, T.J.R. Hughes:
//     Isogeometric finite element data structures based on Bezier extraction of NURBS.
//       International Journal for Numerical Methods in Engineering, Vol. 87(1-5), pp. 15-47.

// Note: indexing was changed to start from 0.

std::vector<MatrixXd> bezierExtractionMatrices(size_t p, const std::vector<double> &knots) {
  std::vector<MatrixXd> C;
  std::vector<double> alphas(p + 1);
  size_t m = knots.size();

  //  Initializations
  size_t a = p + 1, b = a + 1, nb = 0;
  C.push_back(MatrixXd::Identity(p + 1, p + 1));
  while (b < m) {
    C.push_back(MatrixXd::Identity(p + 1, p + 1)); // Initialize the next extraction operator
    size_t i = b;

    // Count multiplicity of the knot at location b
    while (b < m && knots[b] == knots[b-1])
      b++;
    size_t mult = b - i + 1;

    if (mult < p) {
      // Use Eq. (10) in the paper to compute the alphas
      double numer = knots[b-1] - knots[a-1];
      for (size_t j = p; j >= mult + 1; --j)
        alphas[j-mult-1] = numer / (knots[a+j-1] - knots[a-1]);
      size_t r = p - mult;

      // Update the matrix coefficients for r new knots
      for (size_t j = 1; j <= r; ++j) {
        size_t save = r - j;
        size_t s = mult + j;
        for (size_t k = p; k >= s; --k) {
          double alpha = alphas[k-s];
          // The following line corresponds to Eq. (9) in the paper
          C[nb].col(k) = alpha * C[nb].col(k) + (1.0 - alpha) * C[nb].col(k - 1);
        }
        if (b < m) {
          // Update overlapping coefficients of the next operator
          C[nb+1].block(save, save, j + 1, 1) = C[nb].block(p - j, p, j + 1, 1);
        }
      }
      nb++; // Finished with the current operator
      if (b < m) {
        // Update indices for the next operator
        a = b;
        b++;
      }
    }
  }
  return C;
}
