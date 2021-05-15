// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#define NDEBUG 1
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// RjnmfC
arma::field<arma::mat> RjnmfC(arma::mat Xs, arma::mat Xu, int k, double alpha, double lambda, double epsilon, int maxiter, bool verbose);
RcppExport SEXP _conos_RjnmfC(SEXP XsSEXP, SEXP XuSEXP, SEXP kSEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP epsilonSEXP, SEXP maxiterSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xu(XuSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(RjnmfC(Xs, Xu, k, alpha, lambda, epsilon, maxiter, verbose));
    return rcpp_result_gen;
END_RCPP
}
// checkBits
bool checkBits();
RcppExport SEXP _conos_checkBits() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(checkBits());
    return rcpp_result_gen;
END_RCPP
}
// checkOpenMP
bool checkOpenMP();
RcppExport SEXP _conos_checkOpenMP() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(checkOpenMP());
    return rcpp_result_gen;
END_RCPP
}
// cpcaF
Rcpp::List cpcaF(const arma::cube& cov, const arma::vec& ng, int ncomp, int maxit, double tol, Nullable<NumericMatrix> eigenvR, bool verbose);
RcppExport SEXP _conos_cpcaF(SEXP covSEXP, SEXP ngSEXP, SEXP ncompSEXP, SEXP maxitSEXP, SEXP tolSEXP, SEXP eigenvRSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type cov(covSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ng(ngSEXP);
    Rcpp::traits::input_parameter< int >::type ncomp(ncompSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericMatrix> >::type eigenvR(eigenvRSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(cpcaF(cov, ng, ncomp, maxit, tol, eigenvR, verbose));
    return rcpp_result_gen;
END_RCPP
}
// greedyModularityCutC
Rcpp::List greedyModularityCutC(arma::imat& merges, arma::vec& deltaM, int N, int minsize, arma::ivec& labels, double minbreadth, bool flatCut);
RcppExport SEXP _conos_greedyModularityCutC(SEXP mergesSEXP, SEXP deltaMSEXP, SEXP NSEXP, SEXP minsizeSEXP, SEXP labelsSEXP, SEXP minbreadthSEXP, SEXP flatCutSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type merges(mergesSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type deltaM(deltaMSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type minsize(minsizeSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< double >::type minbreadth(minbreadthSEXP);
    Rcpp::traits::input_parameter< bool >::type flatCut(flatCutSEXP);
    rcpp_result_gen = Rcpp::wrap(greedyModularityCutC(merges, deltaM, N, minsize, labels, minbreadth, flatCut));
    return rcpp_result_gen;
END_RCPP
}
// treeJaccard
Rcpp::List treeJaccard(arma::imat& merges, arma::imat& clusters, arma::ivec& clusterTotals, Rcpp::Nullable<arma::imat&> clmerges);
RcppExport SEXP _conos_treeJaccard(SEXP mergesSEXP, SEXP clustersSEXP, SEXP clusterTotalsSEXP, SEXP clmergesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type merges(mergesSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type clusters(clustersSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type clusterTotals(clusterTotalsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::imat&> >::type clmerges(clmergesSEXP);
    rcpp_result_gen = Rcpp::wrap(treeJaccard(merges, clusters, clusterTotals, clmerges));
    return rcpp_result_gen;
END_RCPP
}
// scoreTreeConsistency
Rcpp::List scoreTreeConsistency(arma::imat& test, arma::imat& ref, arma::ivec& leafidmap, int minsize);
RcppExport SEXP _conos_scoreTreeConsistency(SEXP testSEXP, SEXP refSEXP, SEXP leafidmapSEXP, SEXP minsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type test(testSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type ref(refSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type leafidmap(leafidmapSEXP);
    Rcpp::traits::input_parameter< int >::type minsize(minsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(scoreTreeConsistency(test, ref, leafidmap, minsize));
    return rcpp_result_gen;
END_RCPP
}
// maxStableClusters
Rcpp::List maxStableClusters(arma::imat& merges, arma::vec& thresholds, double minthreshold, int minsize);
RcppExport SEXP _conos_maxStableClusters(SEXP mergesSEXP, SEXP thresholdsSEXP, SEXP minthresholdSEXP, SEXP minsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type merges(mergesSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< double >::type minthreshold(minthresholdSEXP);
    Rcpp::traits::input_parameter< int >::type minsize(minsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(maxStableClusters(merges, thresholds, minthreshold, minsize));
    return rcpp_result_gen;
END_RCPP
}
// pareDownHubEdges
NumericVector pareDownHubEdges(SEXP sY, IntegerVector rowN, int k, int klow);
RcppExport SEXP _conos_pareDownHubEdges(SEXP sYSEXP, SEXP rowNSEXP, SEXP kSEXP, SEXP klowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sY(sYSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rowN(rowNSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type klow(klowSEXP);
    rcpp_result_gen = Rcpp::wrap(pareDownHubEdges(sY, rowN, k, klow));
    return rcpp_result_gen;
END_RCPP
}
// getSumWeightMatrix
Rcpp::NumericMatrix getSumWeightMatrix(const std::vector<double>& weights, const std::vector<int>& row_inds, const std::vector<int>& col_inds, const std::vector<int>& factor_levels, bool normalize);
RcppExport SEXP _conos_getSumWeightMatrix(SEXP weightsSEXP, SEXP row_indsSEXP, SEXP col_indsSEXP, SEXP factor_levelsSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type row_inds(row_indsSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type col_inds(col_indsSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type factor_levels(factor_levelsSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(getSumWeightMatrix(weights, row_inds, col_inds, factor_levels, normalize));
    return rcpp_result_gen;
END_RCPP
}
// adjustWeightsByCellBalancingC
std::vector<double> adjustWeightsByCellBalancingC(std::vector<double> weights, const std::vector<int>& row_inds, const std::vector<int>& col_inds, const std::vector<int>& factor_levels, Rcpp::NumericMatrix dividers);
RcppExport SEXP _conos_adjustWeightsByCellBalancingC(SEXP weightsSEXP, SEXP row_indsSEXP, SEXP col_indsSEXP, SEXP factor_levelsSEXP, SEXP dividersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type row_inds(row_indsSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type col_inds(col_indsSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type factor_levels(factor_levelsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dividers(dividersSEXP);
    rcpp_result_gen = Rcpp::wrap(adjustWeightsByCellBalancingC(weights, row_inds, col_inds, factor_levels, dividers));
    return rcpp_result_gen;
END_RCPP
}
// referenceWij
Eigen::SparseMatrix<double> referenceWij(const arma::ivec& i, const arma::ivec& j, arma::vec& d, Rcpp::Nullable<Rcpp::IntegerVector> threads, double perplexity);
RcppExport SEXP _conos_referenceWij(SEXP iSEXP, SEXP jSEXP, SEXP dSEXP, SEXP threadsSEXP, SEXP perplexitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type j(jSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type d(dSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< double >::type perplexity(perplexitySEXP);
    rcpp_result_gen = Rcpp::wrap(referenceWij(i, j, d, threads, perplexity));
    return rcpp_result_gen;
END_RCPP
}
// as_factor
Rcpp::List as_factor(const std::vector<std::string>& vals);
RcppExport SEXP _conos_as_factor(SEXP valsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type vals(valsSEXP);
    rcpp_result_gen = Rcpp::wrap(as_factor(vals));
    return rcpp_result_gen;
END_RCPP
}
// get_nearest_neighbors
Rcpp::List get_nearest_neighbors(const std::vector<std::vector<int>>& adjacency_list, const std::vector<std::vector<double>>& transition_probabilities, int n_verts, int n_cores, double min_prob, int min_visited_verts, double min_prob_lower, int max_hitting_nn_num, int max_commute_nn_num, bool verbose);
RcppExport SEXP _conos_get_nearest_neighbors(SEXP adjacency_listSEXP, SEXP transition_probabilitiesSEXP, SEXP n_vertsSEXP, SEXP n_coresSEXP, SEXP min_probSEXP, SEXP min_visited_vertsSEXP, SEXP min_prob_lowerSEXP, SEXP max_hitting_nn_numSEXP, SEXP max_commute_nn_numSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::vector<int>>& >::type adjacency_list(adjacency_listSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<double>>& >::type transition_probabilities(transition_probabilitiesSEXP);
    Rcpp::traits::input_parameter< int >::type n_verts(n_vertsSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    Rcpp::traits::input_parameter< double >::type min_prob(min_probSEXP);
    Rcpp::traits::input_parameter< int >::type min_visited_verts(min_visited_vertsSEXP);
    Rcpp::traits::input_parameter< double >::type min_prob_lower(min_prob_lowerSEXP);
    Rcpp::traits::input_parameter< int >::type max_hitting_nn_num(max_hitting_nn_numSEXP);
    Rcpp::traits::input_parameter< int >::type max_commute_nn_num(max_commute_nn_numSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(get_nearest_neighbors(adjacency_list, transition_probabilities, n_verts, n_cores, min_prob, min_visited_verts, min_prob_lower, max_hitting_nn_num, max_commute_nn_num, verbose));
    return rcpp_result_gen;
END_RCPP
}
// sgd
arma::mat sgd(arma::mat& coords, arma::ivec& targets_i, arma::ivec& sources_j, arma::ivec& ps, arma::vec& weights, const double& gamma, const double& rho, const arma::uword& n_samples, const int& M, const double& alpha, const Rcpp::Nullable<Rcpp::NumericVector> momentum, const bool& useDegree, const Rcpp::Nullable<Rcpp::NumericVector> seed, const Rcpp::Nullable<Rcpp::IntegerVector> threads, const bool verbose);
RcppExport SEXP _conos_sgd(SEXP coordsSEXP, SEXP targets_iSEXP, SEXP sources_jSEXP, SEXP psSEXP, SEXP weightsSEXP, SEXP gammaSEXP, SEXP rhoSEXP, SEXP n_samplesSEXP, SEXP MSEXP, SEXP alphaSEXP, SEXP momentumSEXP, SEXP useDegreeSEXP, SEXP seedSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type targets_i(targets_iSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type sources_j(sources_jSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type ps(psSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericVector> >::type momentum(momentumSEXP);
    Rcpp::traits::input_parameter< const bool& >::type useDegree(useDegreeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericVector> >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::IntegerVector> >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(sgd(coords, targets_i, sources_j, ps, weights, gamma, rho, n_samples, M, alpha, momentum, useDegree, seed, threads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// propagate_labels
Rcpp::NumericMatrix propagate_labels(const Rcpp::StringMatrix& edge_verts, const std::vector<double>& edge_weights, const Rcpp::StringVector& vert_labels, int max_n_iters, bool verbose, double diffusion_fading, double diffusion_fading_const, double tol, bool fixed_initial_labels);
RcppExport SEXP _conos_propagate_labels(SEXP edge_vertsSEXP, SEXP edge_weightsSEXP, SEXP vert_labelsSEXP, SEXP max_n_itersSEXP, SEXP verboseSEXP, SEXP diffusion_fadingSEXP, SEXP diffusion_fading_constSEXP, SEXP tolSEXP, SEXP fixed_initial_labelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringMatrix& >::type edge_verts(edge_vertsSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type edge_weights(edge_weightsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type vert_labels(vert_labelsSEXP);
    Rcpp::traits::input_parameter< int >::type max_n_iters(max_n_itersSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type diffusion_fading(diffusion_fadingSEXP);
    Rcpp::traits::input_parameter< double >::type diffusion_fading_const(diffusion_fading_constSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type fixed_initial_labels(fixed_initial_labelsSEXP);
    rcpp_result_gen = Rcpp::wrap(propagate_labels(edge_verts, edge_weights, vert_labels, max_n_iters, verbose, diffusion_fading, diffusion_fading_const, tol, fixed_initial_labels));
    return rcpp_result_gen;
END_RCPP
}
// smooth_count_matrix
SEXP smooth_count_matrix(const Rcpp::StringMatrix& edge_verts, const std::vector<double>& edge_weights, const Rcpp::NumericMatrix& count_matrix, const std::vector<bool>& is_label_fixed, int max_n_iters, double diffusion_fading, double diffusion_fading_const, double tol, bool verbose, bool normalize);
RcppExport SEXP _conos_smooth_count_matrix(SEXP edge_vertsSEXP, SEXP edge_weightsSEXP, SEXP count_matrixSEXP, SEXP is_label_fixedSEXP, SEXP max_n_itersSEXP, SEXP diffusion_fadingSEXP, SEXP diffusion_fading_constSEXP, SEXP tolSEXP, SEXP verboseSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringMatrix& >::type edge_verts(edge_vertsSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type edge_weights(edge_weightsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type count_matrix(count_matrixSEXP);
    Rcpp::traits::input_parameter< const std::vector<bool>& >::type is_label_fixed(is_label_fixedSEXP);
    Rcpp::traits::input_parameter< int >::type max_n_iters(max_n_itersSEXP);
    Rcpp::traits::input_parameter< double >::type diffusion_fading(diffusion_fadingSEXP);
    Rcpp::traits::input_parameter< double >::type diffusion_fading_const(diffusion_fading_constSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(smooth_count_matrix(edge_verts, edge_weights, count_matrix, is_label_fixed, max_n_iters, diffusion_fading, diffusion_fading_const, tol, verbose, normalize));
    return rcpp_result_gen;
END_RCPP
}
// adjacent_vertices
Rcpp::List adjacent_vertices(const Rcpp::StringMatrix& edge_verts);
RcppExport SEXP _conos_adjacent_vertices(SEXP edge_vertsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringMatrix& >::type edge_verts(edge_vertsSEXP);
    rcpp_result_gen = Rcpp::wrap(adjacent_vertices(edge_verts));
    return rcpp_result_gen;
END_RCPP
}
// adjacent_vertex_weights
Rcpp::List adjacent_vertex_weights(const Rcpp::StringMatrix& edge_verts, const std::vector<double>& edge_weights);
RcppExport SEXP _conos_adjacent_vertex_weights(SEXP edge_vertsSEXP, SEXP edge_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringMatrix& >::type edge_verts(edge_vertsSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type edge_weights(edge_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(adjacent_vertex_weights(edge_verts, edge_weights));
    return rcpp_result_gen;
END_RCPP
}
// spcov
Eigen::MatrixXd spcov(const Eigen::SparseMatrix<double>& m, Eigen::VectorXd cm);
RcppExport SEXP _conos_spcov(SEXP mSEXP, SEXP cmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type m(mSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type cm(cmSEXP);
    rcpp_result_gen = Rcpp::wrap(spcov(m, cm));
    return rcpp_result_gen;
END_RCPP
}

RcppExport void adjustedRand(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_conos_RjnmfC", (DL_FUNC) &_conos_RjnmfC, 8},
    {"_conos_checkBits", (DL_FUNC) &_conos_checkBits, 0},
    {"_conos_checkOpenMP", (DL_FUNC) &_conos_checkOpenMP, 0},
    {"_conos_cpcaF", (DL_FUNC) &_conos_cpcaF, 7},
    {"_conos_greedyModularityCutC", (DL_FUNC) &_conos_greedyModularityCutC, 7},
    {"_conos_treeJaccard", (DL_FUNC) &_conos_treeJaccard, 4},
    {"_conos_scoreTreeConsistency", (DL_FUNC) &_conos_scoreTreeConsistency, 4},
    {"_conos_maxStableClusters", (DL_FUNC) &_conos_maxStableClusters, 4},
    {"_conos_pareDownHubEdges", (DL_FUNC) &_conos_pareDownHubEdges, 4},
    {"_conos_getSumWeightMatrix", (DL_FUNC) &_conos_getSumWeightMatrix, 5},
    {"_conos_adjustWeightsByCellBalancingC", (DL_FUNC) &_conos_adjustWeightsByCellBalancingC, 5},
    {"_conos_referenceWij", (DL_FUNC) &_conos_referenceWij, 5},
    {"_conos_as_factor", (DL_FUNC) &_conos_as_factor, 1},
    {"_conos_get_nearest_neighbors", (DL_FUNC) &_conos_get_nearest_neighbors, 10},
    {"_conos_sgd", (DL_FUNC) &_conos_sgd, 15},
    {"_conos_propagate_labels", (DL_FUNC) &_conos_propagate_labels, 9},
    {"_conos_smooth_count_matrix", (DL_FUNC) &_conos_smooth_count_matrix, 10},
    {"_conos_adjacent_vertices", (DL_FUNC) &_conos_adjacent_vertices, 1},
    {"_conos_adjacent_vertex_weights", (DL_FUNC) &_conos_adjacent_vertex_weights, 2},
    {"_conos_spcov", (DL_FUNC) &_conos_spcov, 2},
    {"adjustedRand", (DL_FUNC) &adjustedRand, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_conos(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
