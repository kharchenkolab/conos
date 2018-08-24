// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

// merges - merge matrix from walktrap -1, clusters - cluster assignment per cell (0:n.clusters-1), clusterTotals - number of cells per cluster
// [[Rcpp::export]]
Rcpp::List findBestClusterThreshold(arma::imat& merges, arma::ivec& clusters, arma::ivec& clusterTotals) {
  int nleafs=merges.n_rows; // number of leafs
  int nclusters=clusterTotals.n_elem;
  
  // allocate answer matrices
  arma::imat tpn(nclusters,nleafs,arma::fill::zeros); // true positive counts
  arma::imat tnn(nclusters,nleafs,arma::fill::zeros); // true negative counts
  arma::mat ot(nclusters,nleafs); // optimal specificity/sensitivity threshold
  arma::imat oti(nclusters,nleafs); // optimal threshold node id
  
  // temp matrices
  arma::ivec onecol=ones<ivec>(nclusters);
  arma::mat tot2(nclusters,2); // temp ot-like matrix for calculating optimal thresholds per node
  
  // iterate through the merge matrix
  for(int i=0;i<merges.n_rows;i++) {
    arma::mat tot(nclusters,3,arma::fill::zeros); // temp ot-like matrix for comparing optimal threshold between nodes
    arma::imat toti(nclusters,3,arma::fill::zeros); // temp oti-like matrix 
    for(int j=0;j<2;j++) {
      int ni= merges(i,j)-nleafs-1;
      // update count matrices
      if(ni<0) { // leaf node
        // adjust tnn: add ones to all entries except for the correct class
        tnn.col(i)+=onecol;
        int ci=clusters[merges(i,j)];
        if(ci>=0) { // not an NA value
          tnn(ci,i)--;
          tpn(ci,i)++;
          tot(ci,j)=1.0/clusterTotals[ci]; // limited by sensitivity (specificity=1)
          toti(ci,j)=merges(i,j);
        }
      } else { // internal node - simply transfer values
        tnn.col(i)+=tnn.col(ni);
        tpn.col(i)+=tpn.col(ni);
        tot.col(j)=ot.col(ni);
        toti.col(j)=oti.col(ni);
      }
    }
    // recalculate specificity/sensitivity for the merged nodes
    tot2.col(0)=conv_to<vec>::from(tpn.col(i))/clusterTotals; // sensitivity
    tot2.col(1)=1.0-conv_to<vec>::from(tnn.col(i))/clusterTotals; // specificity
    tot.col(2)=min(tot2,1); // report minimal of the sensitivity or specificity achieved
    toti.col(2).fill(i);
    uvec mi=index_max(tot,1); // report maximum threshold achieved between leafs and internal nodes
    // record threshold value and index according to the best threshold
    for(int j=0;j<mi.n_elem;j++) { 
      ot(j,i)=tot(j,mi[j]);
      oti(j,i)=toti(j,mi[j]);
    }
  }
  return List::create(Named("threshold")=ot.col(merges.n_rows-1), Named("node")=oti.col(merges.n_rows-1));
  //Named("tpn")=tpn,Named("tnn")=tnn,Named("ot")=ot,Named("oti")=oti);
}



