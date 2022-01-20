
#include <Rcpp.h>

using namespace Rcpp;

// From clues
// https://github.com/cran/clues/blob/master/src/adjustedRand.c

/***************
 Calculate indices to measure the agreement between two partitions. 
 This c program is written by Weiliang Qiu (Dec. 2001).
 Modified by Weiliang Qiu (Feb. 2006)
 
 If you have any problem or comments on this program, please contact Weiliang Qiu via e-mail
 stwxq@channing.harvard.edu
 Many Thanks!
 Feel free to use this c program for your research purpose. 
@article{Milligan:1986,
  author="{Milligan, G. W.} and {Cooper, M. C.}",
  title={A study of the comparability of external criteria for hierarchical
cluster analysis},
  journal={Multivariate Behavioral Research},
  volume={21},
  pages={441-458},
  year={1986}
}
***************/

/***************
  cl1 --- partition 1 of the data set. 'cl1' is a 'n' by 1 vector.
  cl1u --- unique values of elements in 'cl1'. 'cl1u' is a 'm1' by 1 vector.
  cl2 --- partition 2 of the data set. 'cl2' is a 'n' by 1 vector.
  cl2u --- unique values of elements in 'cl2'. 'cl2u' is a 'm2' by 1 vector.
  m1 --- number of clusters in partition 1
  m2 --- number of clusters in partition 2 (m2 can be not equal to m1)
  nn --- number of data points
  fflag = 1 --- Rand index
  fflag = 2 --- Hubert and Arabie's adjusted Rand index
  fflag = 3 --- Morey and Agresti's adjusted Rand index
  fflag = 4 --- Fowlkes and Mallows's index
  fflag = 5 --- Jaccard index
***************/


// [[Rcpp::export]]
double adjustedRandcpp(Rcpp::NumericVector cl1, Rcpp::NumericVector cl1u, Rcpp::NumericVector cl2, Rcpp::NumericVector cl2u, 
    int mm1, int mm2, int nn, int fflag)
{
    
    // formerly 'double *r12'
    // see the R code how this worked in clues, i.e. 
    // multiple calls used to populate vector
    double result = 0.0;

    int i, j, t, r, *nmatrix;
    //int mm1, mm2, nn, fflag;
    double a, b, c, d, denom; 
    double *nc, *nr, ni_2, n_j2, nt, n_c;
 
    //mm1 = *m1; mm2 = *m2; nn = *n; fflag = *flag;
 
    nmatrix = (int *)malloc((size_t)(mm1 * mm2 * sizeof(int)));
    nc = (double *)malloc((size_t)(mm2 * sizeof(double)));
    nr = (double *)malloc((size_t)(mm1 * sizeof(double)));
 
    a = 0.0; b = 0.0; c = 0.0; d = 0.0;
    denom = 0.0;
    for(t = 0; t < nn ; t ++){
        for(r = t + 1; r < nn; r ++){
            if((cl1[t] == cl1[r]) && (cl2[t] == cl2[r])){
                a = a + 1.0;
            } else if((cl1[t] == cl1[r]) && (cl2[t] != cl2[r])){
                b = b + 1.0;
            } else if((cl1[t] != cl1[r]) && (cl2[t] == cl2[r])){
                c = c + 1.0;
            } else{
                d = d + 1.0;
            }
        }
    }
    // get nij
    for(t = 0; t < mm1; t ++){
        for(r = 0; r < mm2; r ++){
            nmatrix[t * mm2 + r] = 0;
            for(i = 0; i < nn; i ++){
                if((cl1[i] == cl1u[t]) && (cl2[i] == cl2u[r])){
                    nmatrix[t * mm2 + r] += 1;
                }
            }
        }
    }
 
    /* nij_2= \sum_{i=1}^{m_1}\sum_{j=1}^{m_2} n_{ij}^2 */
    /* nr[i]= \sum_{j=1}^{m_2} n_{ij} */
    double nij_2 = 0.0;
    for(i = 0; i < mm1; i ++){
        nr[i] = 0; 
        for(j = 0; j < mm2; j ++){
            nr[i] += nmatrix[i * mm2 + j];
            nij_2 += pow(nmatrix[i * mm2 + j],2);
        }
    }
 
    /* nc[j]= \sum_{i=1}^{m_1} n_{ij} */
    for(i = 0; i < mm2; i ++){
        nc[i] = 0; 
        for(j = 0; j < mm1; j ++){
            nc[i] += nmatrix[j * mm2 + i];
        }
    }
 
    /* ni_2=\sum_{i=1}^{m_1} n_{i.}^2 */
    /* nt=\sum_{i=1}^{m_1}\sum_{j=1}^{m_2} n_{ij} */
    ni_2 = 0.0; n_j2 = 0.0; nt = 0.0;
    for(i = 0; i < mm1; i ++){
        nt += nr[i];
        ni_2 += nr[i] * nr[i];
    }
    /* n_j2=\sum_{j=1}^{m_2} n_{.j}^2 */
    for(j = 0; j < mm2; j ++){
        n_j2 += nc[j] * nc[j];
    }
 
    if(fflag == 2){ //Hubert and Arabie
        n_c = ( nt * (nt * nt + 1.0) - (nt + 1.0) * (ni_2 + n_j2) + 2.0 * ni_2 * n_j2 / nt ) / (2.0 * (nt - 1.0));
        //numer = a + d - n_c;
        denom = a + b + c + d - n_c;
        if(denom < 1.0e-10){ 
            result = 1.0; 
        } else { 
            result = (a + d - n_c) / (a + b + c + d - n_c); 
        }
    } else if(fflag == 3) { //Morey and Agresti
        n_c = nt * (nt - 1.0) / 2.0 - (ni_2 + n_j2) / 2.0 + ni_2 * n_j2 / (nt * nt);
        // numer = a + d - n_c;
        denom = a + b + c + d - n_c;
        if(denom < 1.0e-10){ 
            result = 1.0; 
        } else { 
            result = (a + d - n_c) / (a + b + c + d - n_c); 
        }
    } else if(fflag == 1) { // Rand
        result = (a + d) / (a + b + c + d);
    } else if(fflag == 4) { //Fowlkes and Mallows
        result = a / sqrt((a + b) * (a + c));
    } else if(fflag == 5) { //Jaccard
        result = a / (a + b + c);
    }

    free(nmatrix);
    free(nc);
    free(nr);

    return result;


}