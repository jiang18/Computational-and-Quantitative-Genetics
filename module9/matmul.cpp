#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// C++ function for matrix multiplication using Eigen
// [[Rcpp::export]]
Eigen::MatrixXd matmult(const Eigen::Map<Eigen::MatrixXd>& A, const Eigen::Map<Eigen::MatrixXd>& B, int n_cores) {
  omp_set_num_threads(n_cores);
  return A * B;
}
