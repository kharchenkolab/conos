// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <stack>
#include <queue>

//#define DEBUG 1
#undef DEBUG

using namespace std;
using namespace arma;
using namespace Rcpp;

// get a list of leafs for a given node non-recursively
// returns a vector of leaf indices
vector<int> get_leafs_nr(arma::imat& merges,int node) {
  vector<int> v;
  int nleafs=merges.n_rows+1; // number of leafs in a tree
  //cout<<"nleafs="<<nleafs<<endl;
  //cout<<"g("<<node<<"=="<<(node-nleafs)<<"):"<<flush;
  v.reserve(max(1,node-nleafs+3));
  stack<int> s; s.push(node);
  while(!s.empty()) {
    int j=s.top(); s.pop();
    //cout<<"j="<<j<<' '<<flush;
    if(j<nleafs) { // leaf node
      v.push_back(j);
    } else {
      s.push(merges(j-nleafs,0)); s.push(merges(j-nleafs,1));
    }
  }
  //cout<<endl;
  return(v);
}

typedef std::unordered_map<int, double> IDM;
typedef std::unordered_map<int, arma::ivec > AIVM;
// construct a cluster membership bool vectors for each node participating in the top N splits of the tree in the results hash
// note: currentSplit counts the index from the bottom of merges (i.e. number of splits from the top)
arma::ivec get_top_split_membership(arma::imat& merges, AIVM& results, int N, int currentSplit=0, bool oneMore=false) {
  int nleafs=merges.n_rows+1; // number of vertices is equal to the number of leafs in the tree
  arma::ivec l(nleafs,arma::fill::zeros);
  int currentRow=merges.n_rows-currentSplit-1;

  //cout<<currentSplit<<' '<<(currentRow+nleafs); if(oneMore) { cout<<"+"; }; cout<<" "<<flush<<endl;
  if(currentRow<0) { throw(Rcpp::exception("currentSplit is out of range",currentSplit)); }

  if(currentSplit>=N && !oneMore) { // we're at the last level, use non-recursive function to get the memberships
    vector<int> v=get_leafs_nr(merges,currentRow+nleafs);
    for (auto i : v) { l[i]=true; } // mark cluster vertices
  } else { // use a recursive call, combine the results
    if(merges(currentRow,0)<nleafs) { // leaf node
      l[merges(currentRow,0)]=true;
      results[merges(currentRow,0)]=l;
      if(merges(currentRow,1)<nleafs) { // two leaf nodes
	l[merges(currentRow,1)]=true;
	arma::ivec l2(nleafs,arma::fill::zeros);
	l2[merges(currentRow,1)]=true;
	results[merges(currentRow,1)]=l2;
      } else {
	arma::ivec l2=get_top_split_membership(merges,results,N,merges.n_rows-merges(currentRow,1)+nleafs-1,oneMore=(currentSplit<N));
	l=sign(l+l2);
      }
    } else { // child1 is an internal node
      arma::ivec l2=get_top_split_membership(merges,results,N,merges.n_rows-merges(currentRow,0)+nleafs-1,oneMore=(currentSplit<N));
      if(merges(currentRow,1)<nleafs) { // child2 is a leaf node
	l[merges(currentRow,1)]=true;
	results[merges(currentRow,1)]=l;
	// fold in
	l=sign(l+l2);
      } else { // two internal children
	arma::ivec l3=get_top_split_membership(merges,results,N,merges.n_rows-merges(currentRow,1)+nleafs-1,oneMore=(currentSplit<N));
	l=sign(l2+l3);
      }
    }
  }
  results[currentRow+nleafs]=l;
  return(l);
}

// a quick helper function to fetch cells from vlabs, looking up membership of missing nodes if needed
arma::ivec get_membership(int node, arma::imat& merges, AIVM& vlabs) {
  arma::ivec labs;
  int nleafs=merges.n_rows+1; // number of leafs in the tree
  auto labit=vlabs.find(node);
  if(labit==vlabs.end()) { // construct
    //cout<<"fetching memberships for "<<node<<endl;
    arma::ivec l(nleafs,arma::fill::zeros);
    vector<int> v=get_leafs_nr(merges,node);
    for (auto i : v) { l[i]=true; } // mark cluster vertices
    labs=vlabs[ node ]=l;
  } else {
    labs=labit->second;
  }
  return(labs);
}

// calculate normalized entropy scale for a given set of labels
double normalized_entropy(arma::ivec& labels,int nlabels) {
  // count frequencies
  arma::ivec counts(nlabels,arma::fill::zeros);
  for(auto i=labels.begin(); i!=labels.end(); ++i) {
    counts[ *i ]++;
  }
  double ent=0.0;
  for(auto i=counts.begin(); i!=counts.end(); ++i) {
    if(*i >0) {
      if(*i == labels.n_elem) { return(0.0); } // extreme case
      double p=((double)*i)/((double)labels.n_elem);
      ent-=p * log(p);
    }
  }
  return(ent/log((double)nlabels));
}

// a similar helper function to find node breadth in cache, or calculate it
double get_breadth(int node, arma::ivec& labels, int nlabels, IDM& breadth_map, arma::imat& merges, AIVM& vlabs) {
  auto bit=breadth_map.find(node);
  double breadth;
  if(bit==breadth_map.end()) { // missing, calculate
    arma::ivec membership=get_membership(node,merges,vlabs);
    arma::ivec cllabs=labels.elem(find(membership));
    breadth_map[node]=breadth=normalized_entropy(cllabs,nlabels);
  } else {
    breadth=bit->second;
  }
  return(breadth);
}




// [[Rcpp::export]]
Rcpp::List greedyModularityCutC(arma::imat& merges, arma::vec& deltaM, int N, int minsize, arma::ivec& labels, double minbreadth, bool flatCut) {
  int nleafs=merges.n_rows+1; // number of leafs in the tree
  // a queue to maintain available cuts
  auto comp =[](const std::pair<int,double> &left, const std::pair<int,double> &right) {
    return left.second < right.second;
  };
  priority_queue<pair<int,double>,vector<pair<int,double>>, decltype(comp) > pq(comp);
  AIVM vlabs; // a place to store calculated (cell) memberships for each node
  int nlabels=max(labels)+1;
  IDM breadth_map; // a map to keep breadth results
  //cout<<"pre-calculating membership vectors "<<flush;
  arma::ivec l=get_top_split_membership(merges,vlabs,N);
  //cout<<"done ("<<vlabs.size()<<" nodes)"<<endl;
  //cout<<"["; for(auto kv : vlabs) { cout<<kv.first<<" "; };   cout<<"]"<<endl;

#ifdef DEBUG
  cout<<"walking down the tree ";
#endif

  unordered_set<int> leafNodes;
  vector<int> splitsequence; splitsequence.reserve(N); // sequence of chosen splits
  vector<double> deltamod; deltamod.reserve(N); // keep track of the modularity changes

  if(flatCut) {
    // simply take first N splits
    int n=merges.n_rows-N-1; if(n<0) { n=0; }
    for(int i=merges.n_rows-1; i>n; i--) {
      splitsequence.push_back(i+nleafs);
      deltamod.push_back(deltaM[i]);
    }
  } else {
    // greedy cut
    // start with the top split
    pq.push(pair<int,double>(merges.n_rows -1, deltaM[merges.n_rows -1]));
    int nsplits=0;
    while(!pq.empty() && nsplits<N) {
      int split=pq.top().first;
      double deltam=pq.top().second;
      pq.pop();
#ifdef DEBUG
      cout<<"considering split "<<split<<" ("<<(split+nleafs-2)<<"): ";
#endif
      // record split
      arma::irowvec r=merges.row(split);
      bool take=true;

      // we need to test the split first, for that we'll need to look at the children
      for(int i=0;i<2;i++) {
	// todo: check other requirements, such as min size, breadth, etc.
	if(minsize>0 || minbreadth>0) {
	  arma::ivec membership=get_membership(r[i],merges,vlabs);
	  if(minsize>0 && sum(membership) < minsize) { // size restriction
#ifdef DEBUG
	    cout<<"size "<<minsize<<" < minsize"<<endl;
#endif
	    take=false;
	  }
	  if(take && minbreadth>0) { // check breadth restriction
	    arma::ivec cllabs=labels.elem(find(membership));
	    double breadth=normalized_entropy(cllabs,nlabels);
	    breadth_map[r[i]]=breadth;
	    if(breadth<minbreadth) {
#ifdef DEBUG
	      cout<<"breadth "<<breadth<<" < minbreadth"<<endl;
#endif
	      take=false;
	    }
	  }
	}
      }

      if(take) {  // take the split
#ifdef DEBUG
	cout<<"accepting"<<endl;
#endif
	splitsequence.push_back(split+nleafs);
	leafNodes.erase(split+nleafs); // splitting the node, so it's not a leaf
	deltamod.push_back(deltam);
	// add children to the queue, without checking anything other the fact that they're not leafs
	for(int i=0;i<2;i++) {
	  leafNodes.insert(r[i]);
	  if(r[i]>=nleafs) {
	    pq.push(pair<int,double>(r[i]-nleafs,deltaM[ r[i]-nleafs ]));
	  }
	}
	nsplits++;
      }
    }
  }

  if(splitsequence.empty()) {
#ifdef DEBUG
    cout<<"returning empty sequence"<<endl;
#endif
    return List::create(Named("splitsequence")=splitsequence);
  }

  //cout<<"done"<<endl;
  leafNodes.clear();

  // construct new merge matrix by walking up the merge sequence
  arma::imat nm(splitsequence.size(),2); // new merge matrix
  int nmleafs=nm.n_rows+1; // number of leafs in the new merge matrix
  arma::imat leafContent(nleafs,nmleafs);
  // initialize id counters
  int leafid=0;
  int internalid=nmleafs;
  arma::vec nbreadth(nm.n_rows+nmleafs);
  unordered_map<int,int> idm; // original node id to new node id translation table
#ifdef DEBUG
  cout<<"constructing new merge matrix "<<flush;
#endif
  for(auto it=splitsequence.rbegin(); it!=splitsequence.rend(); ++it) {
    double breadth=get_breadth(*it,labels,nlabels,breadth_map, merges, vlabs);
    int id=idm[*it]=internalid++;
    nbreadth[id]=breadth;
    //cout<<(*it)<<"->"<<id<<" "<<flush;
    arma::irowvec r=merges.row(*it-nleafs);
    // examine content and translate ids
    for(int i=0;i<2;i++) {
      //cout<<" "<<r[i];
      auto si=idm.find(r[i]);
      if(si !=idm.end()) { // internal node, look up id
        //cout<<(r[i])<<"=="<<(si->second)<<" "<<flush;
        double breadth=get_breadth(r[i],labels,nlabels,breadth_map, merges, vlabs);
        r[i]=si->second;
        nbreadth[r[i]]=breadth;
      } else { // leaf node
        //cout<<(r[i])<<"~>"<<leafid<<" "<<flush;
        double breadth=get_breadth(r[i],labels,nlabels,breadth_map, merges, vlabs);
        leafNodes.insert(r[i]);
        int ni=r[i]=idm[r[i]]=leafid++;
        nbreadth[r[i]]=breadth;
      }
    }
    nm.row(id-nmleafs)=r;
  }
#ifdef DEBUG
  cout<<endl;
#endif

  //cout<<"filling out leaf content "<<flush;
  for(auto ln: leafNodes) {
    leafContent.col(idm[ln])=get_membership(ln,merges,vlabs);
  }

  //cout<<"done"<<endl;
  // cout<<"translating ids"<<flush;
  // // translate ids
  // arma::imat oldnm=nm;
  // nm.transform( [&idm](int val) { return (idm[ val ]); } );
  // cout<<endl;

  return List::create(Named("merges")=nm, Named("leafContent")=leafContent,Named("deltaM")=deltamod,
		      Named("splitsequence")=splitsequence,Named("breadth")=nbreadth);

}


// internal: calculate Jaccard coefficients for the leafs and the merge nodes given
// true positive matrix for leaves only and total number of cell per leaf as supplied, as well as cluster totals
Rcpp::List _treeJaccard(arma::imat& merges, arma::imat& tpm, arma::ivec& ncells, arma::ivec& clusterTotals) {
  int nclusters=tpm.n_rows; int nmerges=merges.n_rows; int nleafs=tpm.n_cols;
  if(nleafs!=nmerges+1) stop("nleafs!=nmerges+1");
  // count matries for the merges
  arma::imat mtpm(nclusters,nmerges,arma::fill::zeros);
  arma::ivec mncells(nmerges,arma::fill::zeros);
  // keeping track of optimal answers
  arma::mat ot(nclusters,nmerges); // optimal Jaccard coefficient
  arma::imat oti(nclusters,nmerges); // optimal threshold node id
  // go through the merges
  // temp matrices
  arma::ivec onecol=ones<ivec>(nclusters);
  arma::mat tot2(nclusters,2); // temp ot-like matrix for calculating optimal thresholds per node

  // iterate through the merge matrix
  for(int i=0;i<merges.n_rows;i++) {
    arma::mat tot(nclusters,3,arma::fill::zeros); // temp ot-like matrix for comparing optimal threshold between nodes
    arma::imat toti(nclusters,3,arma::fill::zeros); // temp oti-like matrix
    for(int j=0;j<2;j++) {
      int ni= merges(i,j)-nmerges-1;
      // update count matrices
      if(ni<0) { // leaf node
        mtpm.col(i)+=tpm.col(merges(i,j));
        mncells(i)+=ncells(merges(i,j));
        // calculate tot
        tot.col(j)=conv_to<vec>::from(tpm.col(merges(i,j)))/conv_to<vec>::from(clusterTotals + ncells(merges(i,j))-tpm.col(merges(i,j))); // Jaccard coefficient
        toti.col(j).fill( merges(i,j) );
      } else { // internal node - transfer values for downstream comparison
        mtpm.col(i)+=mtpm.col(ni);
        mncells(i)+=mncells(ni);
        tot.col(j)=ot.col(ni);
        toti.col(j)=oti.col(ni);
      }
    }
    // recalculate Jaccard coefficient for the merged node
    tot.col(2)=conv_to<vec>::from(mtpm.col(i))/conv_to<vec>::from(clusterTotals + mncells(i) - mtpm.col(i)); // Jaccard coefficient
    toti.col(2).fill(i);
    uvec mi=index_max(tot,1); // report maximum threshold achieved between the leafs and the internal nodes
    // record threshold value and index according to the best threshold
    for(int j=0;j<mi.n_elem;j++) {
      ot(j,i)=tot(j,mi[j]);
      oti(j,i)=toti(j,mi[j]);
    }
  }

  return List::create(Named("threshold")=ot.col(merges.n_rows-1), Named("node")=oti.col(merges.n_rows-1));
}

// merges - merge matrix from walktrap -1, clusters - cluster counts per terminal node (1 x nealfs for cell leaves - containing integer cluster ids, or n.clusters x nelafs for cluster leaves - containing cell counts per cluster), clusterTotals - number of cells per cluster
// [[Rcpp::export]]
Rcpp::List treeJaccard(arma::imat& merges, arma::imat& clusters, arma::ivec& clusterTotals, Rcpp::Nullable<arma::imat&> clmerges = R_NilValue) {
  int nleafs=merges.n_rows+1; // number of leafs
  int nclusters=clusterTotals.n_elem;
  // tp and ncells for the leafs
  arma::imat tpm(nclusters,nleafs,arma::fill::zeros); // true positive counts
  arma::ivec ncells(nleafs,arma::fill::zeros); // total number of cells per leaf

  if(clusters.n_rows>1) { // tree with leaves corresponding to clusters
    if(nleafs != clusters.n_cols) stop("number of leafs differs from the number of columns in the clusters matrix");
    tpm=clusters;
    ncells=sum(clusters,0).t();
  } else { // tree with leaves corresponding to individual cells
    for(int i=0;i<nleafs;i++) {
      // adjust tnn: add ones to all entries except for the correct class
      ncells(i)++;
      int ci=clusters(0,i);
      if(ci>=0) { // not an NA value
        tpm(ci,i)++;
      }
    }
  }

  if(clmerges.isNotNull()) {
    arma::imat clm=Rcpp::as<arma::imat>(clmerges);
    // append additional "clusters" corresponding to the merges of earlier clustrers
    arma::imat mtpm(clm.n_rows,nleafs,arma::fill::zeros);
    arma::ivec mclT(clm.n_rows,arma::fill::zeros);
    arma::ivec t1(nleafs);
    for(int i=0;i<clm.n_rows;i++) {
      for(int j=0;j<2;j++) {
        int ni= clm(i,j)-clm.n_rows-1;
        // update count matrices
        if(ni<0) { // leaf node, take counts from fpm/tpm
          mtpm.row(i)+=tpm.row(clm(i,j));
          mclT(i)+=clusterTotals(clm(i,j));
        } else { // internal node, take counts from mtpm/mfpm
          mtpm.row(i)+=mtpm.row(ni);
          mclT(i)+=mclT(ni);
        }
      }
    }
    // join rows
    tpm=join_cols(tpm,mtpm);
    clusterTotals=join_cols(clusterTotals,mclT);
  }

  return _treeJaccard(merges, tpm, ncells, clusterTotals);
}




// CODE FOR TREE vs. TREE SCORING
typedef std::unordered_map<int, vector<int> > IIVM;

// a quick helper function to fetch cells from vlabs, looking up membership of missing nodes if needed
vector<int> get_membership_indices(int node, arma::imat& merges, IIVM& vlabs) {
  vector<int> labs;
  int nleafs=merges.n_rows+1; // number of leafs in the tree
  auto labit=vlabs.find(node);
  if(labit==vlabs.end()) { // construct
    //cout<<"fetching memberships for "<<node<<endl;
    labs=get_leafs_nr(merges,node);
  } else {
    labs=labit->second;
  }
  return(labs);
}

// test - merge matrix from walktrap -1, ref - similar merges matrix for the reference against which to compare
// idmap - a vector of id mappings from ref to test ids (of the same length as the number of leafs in ref; 0-based); those that do not appear should be negative
// [[Rcpp::export]]
Rcpp::List scoreTreeConsistency(arma::imat& test, arma::imat& ref, arma::ivec& leafidmap, int minsize=10) {

  IIVM vlabs; // a place to store calculated (cell) memberships for each node
#ifdef DEBUG
  cout<<"phase I ... "<<flush;
#endif
  // phase I: walk up ref tree to form independent clusters of >=mincells
  int nleafs=ref.n_rows+1; // number of leafs
  arma::ivec leafcount(ref.n_rows,arma::fill::zeros); // a vector to count number of leafs under each merge
  arma::ivec lastmerge(ref.n_rows,arma::fill::zeros); // marking whether the node was the last merge
  for(int i=0;i<ref.n_rows;i++) {
    bool domerge=false;
    bool unmerged_child=false;
    // we will merge
    int lc=0;
    for(int j=0;j<2;j++) {
      int ni= ref(i,j)-nleafs;
      // update count matrices
      if(ni<0) { // leaf node
        domerge=true;
        lc++;
      } else { // inner node
        if(leafcount[ni] < minsize) {
          domerge=true;
        }
        if(lastmerge[ni]==0) { unmerged_child=true;}
        lc+=leafcount[ni];
      }
    }
    leafcount[i]=lc;
    domerge=domerge && !unmerged_child;
    if(domerge) {
      for(int j=0;j<2;j++) {
        int ni= ref(i,j)-nleafs;
        // update count matrices
        if(ni>=0) { // leaf node
          lastmerge[ni]=0;
        }
      }
      lastmerge[i]=1;
    }
  }
#ifdef DEBUG
  cout<<"done."<<endl;
#endif
  arma::imat factor(1,test.n_rows+1);
  arma::vec thresholds(ref.n_rows,arma::fill::zeros); // a vector to keep final thresholds in
  bool hadmerges=true;
  int step=0; // to count interations
  // phase II: iterative score/merge (iterate between part 1/2 below)
  // part 1: score with current ref cut
  while(hadmerges) {
    // current cut
    arma::uvec currentCut=find(lastmerge);
#ifdef DEBUG
    cout<<"step "<<step<<": "<<currentCut.n_elem<<" clusters. constructing factor ... "<<flush;
#endif
    factor.fill(-1); // by default, cells are left unclassified
    // make up a factor
    arma::ivec factorTotals(currentCut.n_elem,arma::fill::zeros);
    for(int i=0;i<currentCut.n_elem;i++) {
      vector<int> v=get_membership_indices(currentCut[i]+nleafs,ref,vlabs);
      //if(i<10) cout<<"i="<<i<<" ("<<v.size()<<"):"<<flush;
      int realsize=0;
      for(auto j: v) {
        if(j>=leafidmap.n_elem) throw std::range_error("requesting out of range node in leafidmap");
        int ni=leafidmap[j];
        if(ni> ((int)test.n_rows)) {
          cout<<"idmap violation at j="<<j<<" ("<<ni<<"); test.n_rows="<<test.n_rows<<' '<<(ni>test.n_rows)<<endl<<flush;
          throw std::range_error("out of range value in leafidmap "+j);
        }
        //if(i<10) cout<<j<<"="<<ni<<' '<<flush;
        //if(ni>=0 && factor[ni]>=0) { cout<<"overriding factor["<<ni<<"]="<<(factor[ni])<<" with "<<i<<endl; }

        if(ni>=0) { factor(0, ni )=i; realsize++; }
      }
      factorTotals[i]=realsize;
      //if(i<10) cout<<endl;
    }
#ifdef DEBUG
    cout<<"done; ";

    cout<<"scoring ... "<<flush;
#endif
    // score current cut
    Rcpp::List a=treeJaccard(test, factor, factorTotals);
    arma::vec currentThresholds = as<arma::vec>(a["threshold"]);


    // persist thresholds
    for(int i=0;i<currentCut.n_elem;i++) {
      thresholds[currentCut[i]]=max(thresholds[currentCut[i]],currentThresholds[i]);
    }



#ifdef DEBUG
    cout<<"merging ... ";
#endif

    // step 2: walk up ref tree to merge all independent clusters (i.e. not hitting previously merged nodes, and not merging with unprocesssed nodes)
    hadmerges=false;
    int nmerges=0;
    arma::ivec justmerged(ref.n_rows,arma::fill::zeros);
    for(int i=0;i<ref.n_rows;i++) {
      if(lastmerge[i]) continue;
      int cmerged=0;
      int cleafs=0;
      for(int j=0;j<2;j++) {
        int ni= ref(i,j)-nleafs;
        if(ni<0) {
          cleafs++;
        } else {
          if(lastmerge[ni] && !justmerged[ni]) cmerged++;
        }
      }
      // merge if one child is a merged node, and
      //   - the other child is a leaf
      //   - the other child is another merged node (that wasn't just merged)
      if(cmerged==2 || (cmerged==1 && cleafs==1)) {
        // merge
        for(int j=0;j<2;j++) {
          int ni= ref(i,j)-nleafs;
          // update count matrices
          if(ni>=0) { // inner node
            lastmerge[ni]=0;
          }
        }
        lastmerge[i]=1;
        justmerged[i]=1;
        nmerges++;
      }
    }
#ifdef DEBUG
    cout<<nmerges<<" merges"<<endl;
#endif
    if(nmerges>0) hadmerges=true;
    step++;
    if(step>100) break;
  }

  return List::create(Named("thresholds")=thresholds,
                      Named("leafcount")=leafcount);

}

// given per-merge thresholds, determine max number of independent stable clusters
// [[Rcpp::export]]
Rcpp::List maxStableClusters(arma::imat& merges, arma::vec& thresholds, double minthreshold=0.8, int minsize=10) {

  // phase I: propagate max subtree thresholds up
  arma::vec maxthreshold(merges.n_rows,arma::fill::zeros);
  arma::ivec leafcount(merges.n_rows,arma::fill::zeros); // a vector to count number of leafs under each merge
  arma::ivec wheresmax(merges.n_rows,arma::fill::zeros); //0 - self; 1- right; 2-left;
  int nleafs=merges.n_rows+1; // number of leafs
  arma::vec ct(2);
#ifdef DEBUG
  cout<<"propagating maxthreshold ... "<<flush;
#endif
  for(int i=0;i<merges.n_rows;i++) {
    int lc=0;
    for(int j=0;j<2;j++) {
      int ni= merges(i,j)-nleafs;
      if(ni<0) { // leaf node
        lc++;
        ct[j]=0;
      } else { // inner node
        if(leafcount[ni] < minsize) {
          ct[j]=0;
        } else {
          ct[j]=maxthreshold[ni];
        }
        lc+=leafcount[ni];
      }
    }
    leafcount[i]=lc;
    if(lc>=minsize) {
      double mct=max(ct);
      if(thresholds[i]<mct) { // record where the maximum was
        wheresmax[i]=(ct[1]>ct[0])+1;
        maxthreshold[i]=mct;
      } else {
        maxthreshold[i]=thresholds[i];
      }
    } else { // current node is too small
      maxthreshold[i]=0;
    }
#ifdef DEBUG
    //cout<<i<<":["<<(merges(i,0)-nleafs)<<","<<(merges(i,1)-nleafs)<<"] lc="<<lc<<" th="<<thresholds[i]<<" mth="<<maxthreshold[i]<<endl;
#endif
  }
#ifdef DEBUG
  cout<<"done"<<endl;
  cout<<"determining max clusters ... ";
#endif

  // phase II: walk downwards
  arma::uvec x=find(maxthreshold>=minthreshold);
  vector<int> terminalnodes; terminalnodes.reserve(x.n_elem);
  stack<int> s; s.push(merges.n_rows-1 + nleafs); // start with the top merge
  while(!s.empty()) {
    int j=s.top(); s.pop();
    if(j>=nleafs) { // internal node
      // check if branch if we should continue
      bool added=false;
      for(int i=0;i<2;i++) {
        int ni=merges(j-nleafs,i);
        if(ni>=nleafs) { // internal
          if(maxthreshold[ni-nleafs]  >= minthreshold) {
             s.push(ni); // will try to split this node in the future
             added=true;
          }
        }
      }
      if(!added) { // push  as a terminal node
        terminalnodes.push_back(j);
      }
    }
  }
#ifdef DEBUG
  cout<<"done ("<<terminalnodes.size()<<" clusters)"<<endl;
#endif
  terminalnodes.shrink_to_fit();

  return List::create(Named("terminalnodes")=terminalnodes,
                      Named("maxthreshold")=maxthreshold,
                      Named("wheresmax")=wheresmax);
}

