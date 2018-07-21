// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
// [[Rcpp::export]]
Eigen::MatrixXd spcov(const Eigen::SparseMatrix<double>& m,Eigen::VectorXd cm) {
  Eigen::MatrixXd v=m.transpose()*m; 
  v-=(cm* cm.transpose())*((double)m.rows());
  v/=((double) m.rows()-1);
  return(v);
}
