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
arma::ivec get_top_split_membership(arma::imat& merges,AIVM& results,int N,int currentSplit=0,bool oneMore=false) {
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
Rcpp::List greedyModularityCut(arma::imat& merges, arma::vec& deltaM, int N, int minsize, arma::ivec& labels, double minbreadth) {
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
