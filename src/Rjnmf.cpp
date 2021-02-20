#include <RcppArmadillo.h>
#include <stdlib.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace std;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif


double inline tr(mat A, mat B) {
  return accu(A % B);
}


#define MIN_VALUE (1e-10)


// [[Rcpp::export]]
arma::field<arma::mat> RjnmfC(arma::mat Xs, arma::mat Xu, int k, double alpha, double lambda, double epsilon, int maxiter, bool verbose) {

  int n = Xs.n_rows;
  int v1 = Xs.n_cols;
  int v2 = Xu.n_cols;

  double lastObj = 0;

  // Make deterministic
  //arma_rng::set_seed(5134);

  // Initialise Matrices
  arma::mat W(n, k, fill::randu);
  for (mat::iterator i = W.begin(); i != W.end(); ++i) *i = abs(*i);

  arma::mat Hs(k, v1, fill::randu);
  for (mat::iterator i = Hs.begin(); i != Hs.end(); ++i) *i = abs(*i);

  arma::mat Hu(k, v2, fill::randu);
  for (mat::iterator i = Hu.begin(); i != Hu.end(); ++i) *i = abs(*i);

  double beta = 1.0 - alpha;
  double trXstXs = tr(Xs, Xs);
  double trXutXu = tr(Xu, Xu);

  mat Wt = W.t();
  mat WtW = Wt * W;
  mat WtXs = Wt * Xs;
  mat WtXu = Wt * Xu;
  mat WtWHs = WtW * Hs;
  mat WtWHu = WtW * Hu;

  int itNum = 1;
  double delta = 2 * epsilon;

  while ( (delta > epsilon) && (itNum <= maxiter)) {

    // Update Hs
    mat Hs_a = alpha * WtXs;
    mat Hs_b = alpha * WtWHs + lambda * Hs;
    for (mat::iterator i = Hs_b.begin(); i != Hs_b.end(); ++i) {if(*i < MIN_VALUE) *i = MIN_VALUE;};
    Hs = Hs % (Hs_a / Hs_b);

    // Update Hu
    mat Hu_a = beta * WtXu;
    mat Hu_b = beta * WtWHu + lambda * Hu; // todo cap
    for (mat::iterator i = Hu_b.begin(); i != Hu_b.end(); ++i) {if(*i < MIN_VALUE) *i = MIN_VALUE;};
    Hu = Hu % (Hu_a / Hu_b);

    // Update W
    mat W_a = alpha * Xs * Hs.t() + beta * Xu * Hu.t();
    mat W_b = alpha * W  * Hs * Hs.t() + beta * W * Hu * Hu.t() + lambda * W;
    for (mat::iterator i = W_b.begin(); i != W_b.end(); ++i) {if (*i < MIN_VALUE) *i = MIN_VALUE;};
    W = W % (W_a / W_b);

    // Calculate objective function

    Wt = W.t();
    WtW = Wt * W;
    WtXs = Wt * Xs;
    WtXu = Wt * Xu;
    WtWHs = WtW * Hs;
    WtWHu = WtW * Hu;

    double tr1 = alpha  * (trXstXs - 2*tr(Hs, WtXs) + tr(Hs, WtWHs));
    double tr2 = beta   * (trXutXu - 2*tr(Hu, WtXu) + tr(Hu, WtWHu));
    double tr3 = lambda * (trace(WtW) + tr(Hs, Hs) + tr(Hu, Hu));
    double Obj = tr1 + tr2 + tr3;

    if (itNum != 1) {
        delta = abs(Obj - lastObj);
        if (verbose) Rcout << "Iteration " << itNum << " Objective: " << Obj << " Delta " << delta << endl << flush;
    } else {
        if (verbose) Rcout << "Iteration " << itNum << " Objective: " << Obj << endl << flush;
    }

    lastObj = Obj;

    Rcpp::checkUserInterrupt();

    ++itNum;
  }

  arma::field<arma::mat> returnList(3);
  returnList(0) = Hs;
  returnList(1) = Hu;
  returnList(2) = W;

  return returnList;
}

