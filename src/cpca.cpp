// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <math.h>
using namespace std;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List cpcaF(const arma::cube& cov,const arma::vec& ng,int ncomp=10,int maxit=1000,double tol=1e-6,Nullable<NumericMatrix> eigenvR=R_NilValue,bool verbose=true) {
  int k = cov.n_slices;
  int p = cov.n_rows;
  int n = sum(ng);
  
  mat D(k,ncomp,fill::zeros);
  mat CPC(p,ncomp,fill::zeros);
  mat Qw(p,p,fill::eye);

  bool converged=false;
  vec itComp(ncomp);

  // initialization
  mat S(p,p,fill::zeros);
  for(int i=0;i<k;i++) {
    S+=cov.slice(i) * (ng[i]/((double)n));
  }

  vec eigval; mat eigvec;
  if(eigenvR.isNotNull()) {
    eigvec=as<mat>(eigenvR);
  } else {
    eig_sym(eigval, eigvec, S);  
    // largest eigenvalues will be at teh front
    eigvec=fliplr(eigvec);
  }
  
  if(verbose) cout<<ncomp<<" CPCAs [";
  for(int comp=0; comp<ncomp; comp++) {
    //if(verbose) cout<<"comp "<<(comp+1)<<endl;
    if(verbose) cout<<'.'<<flush;
    vec q(eigvec.col(comp));
    vec d(k,fill::zeros);
    for(int i=0;i<k;i++) {
      arma::vec v=q.t() * (cov.slice(i) * q);
      d[i]=v[0];
    }
    //cout<<d<<endl;
    double cost0=0;
    int it=0;
    for(it=0; it<maxit; it++) {
      S.zeros();
      for(int i=0;i<k;i++) {
	S+=cov.slice(i) * (ng[i]/d[i]);
      }
      vec w=S*q;
      
      if(comp!=0) {
	w=Qw*w;
      }
      
      q=w/sqrt(sum(w%w));
      
      for(int i=0;i<k;i++) {
	arma::vec v=q.t() * (cov.slice(i) * q);
	d[i]=v[0];
      }
       
      double cost=sum(log(d) % ng);
      double delta=fabs((cost-cost0)/cost);
      //cout<<"delta["<<it<<"]="<<delta<<endl;
      if(delta<tol) break;
      cost0 = cost;
    }
    itComp[comp]=it;
    D.col(comp)=d;
    CPC.col(comp)=q;
    //arma::mat v=(q * q.t());
    //cout<<"q*t(q)="<<v<<endl;
    Qw-=q * q.t();
  }
  cout<<"]"<<endl;
  
  // return
  List ret;
  ret["D"] = D.t();
  ret["CPC"] = CPC;
  ret["niter"] = itComp;
  return(ret);

}

