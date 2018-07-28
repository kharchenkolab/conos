#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>

#include <RcppEigen.h>
#include <progress.hpp>

using si_map_t = std::unordered_map<std::string, size_t>;
using s_vec_t = std::vector<std::string>;
using SpMat = Eigen::SparseMatrix<double>;
using Mat = Eigen::MatrixXd;

class Edge {
public:
  const size_t v_start;
  const size_t v_end;
  const double weight;
  const double length;

  Edge(size_t v_start, size_t v_end, double weight)
    : v_start(v_start), v_end(v_end), weight(weight), length(1 - weight)
  {}
};

si_map_t order_strings(const s_vec_t &vec) {
  si_map_t nums;
  for (auto const &s : vec) {
    nums.emplace(s, nums.size());
  }

  return nums;
}

si_map_t parse_edges(const Rcpp::StringMatrix &edge_verts, const std::vector<double> &edge_weights, std::vector<Edge> &edges, s_vec_t vertex_names=s_vec_t()) {
  if (edge_verts.nrow() != edge_weights.size() || edge_verts.ncol() != 2)
    Rcpp::stop("Incorrect dimension of input vectors");

  if (vertex_names.empty()) {
    for (size_t i = 0; i < edge_verts.nrow(); ++i) {
      vertex_names.push_back(Rcpp::as<std::string>(edge_verts(i, 0)));
      vertex_names.push_back(Rcpp::as<std::string>(edge_verts(i, 1)));
    }
  }

  auto vertex_ids = order_strings(vertex_names);
  for (size_t i = 0; i < edge_verts.nrow(); ++i) { // TODO: add informative message in case vertex name is not presented in the map
    edges.emplace_back(vertex_ids.at(Rcpp::as<std::string>(edge_verts(i, 0))), vertex_ids.at(Rcpp::as<std::string>(edge_verts(i, 1))), edge_weights.at(i));
  }

  return vertex_ids;
}

void smooth_count_matrix(const std::vector<Edge> &edges, Mat &count_matrix, int max_n_iters, double diffusion_fading, double diffusion_fading_const) {
  std::vector<double> sum_weights(count_matrix.rows(), 1);
  Progress p(max_n_iters * edges.size(), true);
  double min_weight = 1e10, max_weight = 0;
  for (int iter = 0; iter < max_n_iters; ++iter) {
    Mat cm_new(count_matrix);
    for (auto const &e : edges) {
      if (Progress::check_abort())
        return;

      double weight = exp(-diffusion_fading * (e.length + diffusion_fading_const));
      min_weight = std::min(min_weight, weight);
      max_weight = std::max(max_weight, weight);

      cm_new.row(e.v_start) += count_matrix.row(e.v_end) * weight;
      cm_new.row(e.v_end) += count_matrix.row(e.v_start) * weight;

      sum_weights.at(e.v_start) += weight;
      sum_weights.at(e.v_end) += weight;

      p.increment();
    }

    for (size_t row_id = 0; row_id < cm_new.rows(); ++row_id) {
      cm_new.row(row_id) /= sum_weights.at(row_id);
      sum_weights.at(row_id) = 1;
    }

    double diff = (cm_new - count_matrix).array().abs().matrix().lpNorm<Eigen::Infinity>();
    std::cout << iter << ": " << diff << std::endl;
    count_matrix = cm_new;
  }

  std::cout << "Min weight: " << min_weight << ", max weight: " << max_weight << ", fading: ("
            << diffusion_fading << ", " << diffusion_fading_const << ")" << std::endl;
}

// [[Rcpp::export]]
SEXP smooth_count_matrix(const Rcpp::StringMatrix &edge_verts, const std::vector<double> &edge_weights, const Rcpp::NumericMatrix &count_matrix,
                         int max_n_iters=10, double diffusion_fading=1.0, double diffusion_fading_const=0.1) {
  // Parse data
  Mat cm_eigen(Rcpp::as<Mat>(count_matrix));

  auto const gene_names = Rcpp::as<s_vec_t>(Rcpp::colnames(count_matrix));
  auto const cell_names = Rcpp::as<s_vec_t>(Rcpp::rownames(count_matrix));

  std::vector<Edge> edges;
  auto vertex_ids = parse_edges(edge_verts, edge_weights, edges, cell_names);

  // Process data
  smooth_count_matrix(edges, cm_eigen, max_n_iters, diffusion_fading, diffusion_fading_const);

  // Convert result back to R
  Rcpp::NumericMatrix cm_res = Rcpp::wrap(cm_eigen);
  Rcpp::colnames(cm_res) = Rcpp::colnames(count_matrix);
  Rcpp::rownames(cm_res) = Rcpp::rownames(count_matrix);

  return cm_res;
}

// [[Rcpp::export]]
Rcpp::List adjacent_vertices(const Rcpp::StringMatrix &edge_verts) {
  std::unordered_map<std::string, std::vector<std::string>> adj_verts;
  for (size_t i = 0; i < edge_verts.nrow(); ++i) {
    auto v1 = Rcpp::as<std::string>(edge_verts(i, 0));
    auto v2 = Rcpp::as<std::string>(edge_verts(i, 1));
    adj_verts[v1].push_back(v2);
    adj_verts[v2].push_back(v1);
  }

  return Rcpp::wrap(adj_verts);
}

// [[Rcpp::export]]
Rcpp::List adjacent_vertex_weights(const Rcpp::StringMatrix &edge_verts, const std::vector<double> &edge_weights) {
  std::unordered_map<std::string, std::vector<double>> adj_verts;
  for (size_t i = 0; i < edge_verts.nrow(); ++i) {
    auto v1 = Rcpp::as<std::string>(edge_verts(i, 0));
    auto v2 = Rcpp::as<std::string>(edge_verts(i, 1));
    adj_verts[v1].push_back(edge_weights.at(i));
    adj_verts[v2].push_back(edge_weights.at(i));
  }

  return Rcpp::wrap(adj_verts);
}
