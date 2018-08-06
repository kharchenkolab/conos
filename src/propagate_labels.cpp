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

class Vertex {
public:
  const std::string name;
  std::vector<double> label_weights;

  Vertex(const std::string &name, size_t n_labels)
    : name(name), label_weights(n_labels, 0)
  {}

  std::vector<double> label_probs() const {
    std::vector<double> probs = this->label_weights;
    double weight_sum = std::accumulate(probs.begin(), probs.end(), 0);
    if (weight_sum < 1e-10)
      return probs;

    for (int i = 0; i < probs.size(); ++i) {
      probs[i] /= weight_sum;
    }

    return probs;
  }
};

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

// Label propogation

void propagate_labels(Vertex &v1, Vertex &v2, double weight) {
  auto probs1 = v1.label_probs();
  auto probs2 = v2.label_probs();

  for (size_t l_id = 0; l_id < probs1.size(); l_id ++) {
    v1.label_weights[l_id] += probs2[l_id] * weight;
    v2.label_weights[l_id] += probs1[l_id] * weight;
  }
}

void propagate_labels(const std::vector<Edge> &edges, std::vector<Vertex> &vertices, int max_n_iters, bool verbose) {
  Progress p(max_n_iters, verbose);
  for (int iter = 0; iter < max_n_iters; ++iter) {
    if (Progress::check_abort())
      break;

    for (auto const e : edges) {
      propagate_labels(vertices.at(e.v_start), vertices.at(e.v_end), e.weight);
    }

    p.increment();
  }
}

// [[Rcpp::export]]
Rcpp::NumericMatrix propagate_labels(const Rcpp::StringMatrix &edge_verts, const std::vector<double> &edge_weights, const Rcpp::StringVector &vert_labels, int max_n_iters=10, bool verbose=true) {
  if (edge_verts.nrow() != edge_weights.size() || edge_verts.ncol() != 2)
    Rcpp::stop("Incorrect dimension of input vectors");

  // Parse data
  std::vector<Edge> edges;
  auto const vertex_ids = parse_edges(edge_verts, edge_weights, edges);
  auto const label_ids = order_strings(Rcpp::as<s_vec_t>(vert_labels));

  std::vector<Vertex> vertices;
  for (auto const p : vertex_ids) {
    vertices.emplace_back(p.first, label_ids.size());
  }

  for (size_t i = 0; i < vert_labels.size(); ++i) {
    size_t label_id = label_ids.at(Rcpp::as<std::string>(vert_labels.at(i)));
    size_t vert_id = vertex_ids.at(Rcpp::as<std::string>(Rcpp::as<Rcpp::StringVector>(vert_labels.names()).at(i)));
    vertices.at(vert_id).label_weights.at(label_id) += 1;
  }

  // Process data
  propagate_labels(edges, vertices, max_n_iters, verbose);

  // Convert result back to R
  Rcpp::NumericMatrix res(vertices.size(), label_ids.size());
  for (int v_id = 0; v_id < vertices.size(); ++v_id) {
    for (int l_id = 0; l_id < label_ids.size(); ++l_id) {
      res(v_id, l_id) = vertices.at(v_id).label_weights.at(l_id);
    }

    res(v_id, Rcpp::_) = res(v_id, Rcpp::_) / Rcpp::sum(res(v_id, Rcpp::_));
  }

  s_vec_t labels(label_ids.size());
  for (auto const &p : label_ids) {
    labels.at(p.second) = p.first;
  }

  s_vec_t vert_names(vertex_ids.size());
  for (auto const &p : vertex_ids) {
    vert_names.at(p.second) = p.first;
  }

  Rcpp::colnames(res) = Rcpp::wrap(labels);
  Rcpp::rownames(res) = Rcpp::wrap(vert_names);

  return res;
}

// Smooth count matrix

void smooth_count_matrix(const std::vector<Edge> &edges, Mat &count_matrix, int max_n_iters, double diffusion_fading, double diffusion_fading_const, double tol, bool verbose) {
  std::vector<double> sum_weights(count_matrix.rows(), 1);
  Progress p(max_n_iters * edges.size(), verbose);
  double min_weight = 1e10, max_weight = 0;
  double inf_norm = 1e10;
  int iter = 0;
  for (iter = 0; iter < max_n_iters; ++iter) {
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

    inf_norm = (cm_new - count_matrix).array().abs().matrix().lpNorm<Eigen::Infinity>();
    if (inf_norm < tol)
      break;

    count_matrix = cm_new;
  }

  if (!verbose)
    return;

  std::cout << "Stop after " << iter << " iterations. Norm: " << inf_norm << std::endl
            << "Min weight: " << min_weight << ", max weight: " << max_weight << ", fading: ("
            << diffusion_fading << ", " << diffusion_fading_const << ")" << std::endl;
}

// [[Rcpp::export]]
SEXP smooth_count_matrix(const Rcpp::StringMatrix &edge_verts, const std::vector<double> &edge_weights, const Rcpp::NumericMatrix &count_matrix,
                         int max_n_iters=10, double diffusion_fading=1.0, double diffusion_fading_const=0.1, double tol=1e-3, bool verbose=true) {
  // Parse data
  Mat cm_eigen(Rcpp::as<Mat>(count_matrix));

  auto const cell_names = Rcpp::as<s_vec_t>(Rcpp::rownames(count_matrix));

  std::vector<Edge> edges;
  parse_edges(edge_verts, edge_weights, edges, cell_names);

  // Process data
  smooth_count_matrix(edges, cm_eigen, max_n_iters, diffusion_fading, diffusion_fading_const, tol, verbose);

  // Convert result back to R
  Rcpp::NumericMatrix cm_res = Rcpp::wrap(cm_eigen);
  Rcpp::colnames(cm_res) = Rcpp::colnames(count_matrix);
  Rcpp::rownames(cm_res) = Rcpp::rownames(count_matrix);

  return cm_res;
}

// Functions for debugging

// [[Rcpp::export]]
Rcpp::List adjacent_vertices(const Rcpp::StringMatrix &edge_verts) {
  std::unordered_map<std::string, s_vec_t> adj_verts;
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
