// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;


// calculates factor-stratified sums for each column
// sY - sparse matrix (dgCmatrix)
// rowSel is a proper R factor
// note that the 0-th column will return sums for any NA values; 0 or negative values will be omitted
// [[Rcpp::export]]
NumericMatrix colSumByFactor(SEXP sY,  IntegerVector rowSel) {
  // need to do this as SEXP, modify the slots on the fly
  S4 mat(sY);
  const arma::uvec i(( unsigned int *)INTEGER(mat.slot("i")),LENGTH(mat.slot("i")),false,true);
  const arma::ivec dims(INTEGER(mat.slot("Dim")),LENGTH(mat.slot("Dim")),false,true);
  const arma::ivec p(INTEGER(mat.slot("p")),LENGTH(mat.slot("p")),false,true);
  arma::vec Y(REAL(mat.slot("x")),LENGTH(mat.slot("x")),false,true);

  List dimnames(mat.slot("Dimnames"));
  CharacterVector geneNames(dimnames[1]);
  

  CharacterVector factorLevels=rowSel.attr("levels");
  int nlevels=factorLevels.size();
  CharacterVector expandedFactorLevels(nlevels+1);
  expandedFactorLevels[0]="NA";
  for(int i=0;i<nlevels;i++) { expandedFactorLevels[i+1]=factorLevels[i]; }
  const arma::ivec rs=arma::ivec(INTEGER(rowSel),LENGTH(rowSel),false,true);

  int ncols=p.size()-1;

  if(nlevels==0) { stop("colSumByFactor(): supplied factor doesn't have any levels!"); }
  //arma::mat sumM(nlevels+1,ncols,arma::fill::zeros);
  NumericMatrix sumM(nlevels+1,ncols); 

  // for each gene
  //#pragma omp parallel for shared(sumM) 
  for(int g=0;g<ncols;g++) {
    int p0=p[g]; int p1=p[g+1];
    if(p1-p0 <1) { continue; }
    for(int j=p0;j<p1;j++) {
      int row=i[j];
      int f=rs[row];
      if(f==NA_INTEGER) {
        sumM(0,g)+=Y[j];
      } else if(f>0) {
        sumM(f,g)+=Y[j];
      }
    }
  }

  
  colnames(sumM) = geneNames;
  rownames(sumM) = expandedFactorLevels;
  return sumM;

}

