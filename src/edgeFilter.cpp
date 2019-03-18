// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;


// greedy algorithm for removing edges (sets their @x to 0) from a sparce matrix, starting with higher-conencted nodes
// guarantees not to remove edges from nodes that have <k edges
// sY - sparse matrix (dgCmatrix)
// rowN - number of non-zero elements per row
// k - minimal number of edges to keep
// will return a modified copy of the @x
// [[Rcpp::export]]
NumericVector pareDownHubEdges(SEXP sY,  IntegerVector rowN, int k) {
  S4 mat(sY);
  const arma::uvec i(( unsigned int *)INTEGER(mat.slot("i")),LENGTH(mat.slot("i")),false,true);
  const arma::ivec dims(INTEGER(mat.slot("Dim")),LENGTH(mat.slot("Dim")),false,true);
  const arma::ivec p(INTEGER(mat.slot("p")),LENGTH(mat.slot("p")),false,true);
  NumericVector x(mat.slot("x"));

  int ncols=p.size()-1;
  int j,jl;
  for(int g=0; g<ncols;g++) {
    
    int p0=p[g]; int p1=p[g+1];
    // check how mant edges the node has to begin with
    int nedges=p1-p0;
    if(nedges <= k) { continue; }

    // order edges
    arma::ivec cid(nedges);
    for(j=p0,jl=0; j<p1; j++,jl++) {
      int row=i[j];
      cid[jl]=rowN[row];
    }
    // order by decreasing degree
    arma::uvec cido=sort_index(cid,"descend");
    // remove top edges, adjusting row counts
    j=0;
    while(nedges > k && j<cido.n_elem && cid[ cido[j] ]>k) {
      x[ cido[j]+p0 ]=0;  // zero-out the edge
      rowN[ i[ cido[j]+p0 ] ]--; // decrease row counts
      nedges--; j++; // adjust counters
    }

  }

  return x;

}
