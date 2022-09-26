

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
#define NDEBUG 1
#include <RcppEigen.h>
#include <Rcpp.h>


// [[Rcpp::export]]
Eigen::MatrixXd spcov(const Eigen::SparseMatrix<double>& m,Eigen::VectorXd cm) {
  Eigen::MatrixXd v=m.transpose()*m; 
  v-=(cm* cm.transpose())*((double)m.rows());
  v/=((double) m.rows()-1);
  return(v);
}


// quick matrix correlation function
// [[Rcpp::export]]
arma::mat arma_mat_cor(const arma::mat& m) {
  return(cor(m));
}
