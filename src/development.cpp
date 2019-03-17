#include <Rcpp.h>
#include <progress.hpp>

#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace Rcpp;

class EdgeRemovalInfo {
private:
  int n1;
  int n2;
  double weight;

public:
  int id;

private:
  int min_n_edges() const {
    return std::min(this->n1, this->n2);
  }

  int max_n_edges() const {
    return std::max(this->n1, this->n2);
  }

public:
  EdgeRemovalInfo(int id, int n1, int n2, double weight)
    : id(id), n1(n1), n2(n2), weight(weight)
  {}

  bool operator <(const EdgeRemovalInfo & e2) const {
    int min1 = this->min_n_edges(), min2 = e2.min_n_edges();
    if (min1 == min2) {
      if (std::abs(this->weight - e2.weight) < 0.01)
        return this->max_n_edges() > e2.max_n_edges();

      return this->weight < e2.weight;
    }

    return min1 > min2;
  }
};

//' @export
// [[Rcpp::export]]
std::vector<bool> edge_removal_mask(const std::vector<std::string> &verts1, const std::vector<std::string> &verts2, const std::vector<double> &weights, int min_neighb_per_vertex, bool verbose=true) {
  if (verts1.size() != verts2.size() || verts1.size() != weights.size())
    stop("Vertex vectors must have the same size");

  Progress p(verts1.size(), false);

  if (verbose) {
    std::cout << "Estimating number of edges per vertex..." << std::endl;
  }

  std::unordered_map<std::string, int> n_neighb_map;
  for (size_t i = 0; i < verts1.size(); ++i) {
    n_neighb_map[verts1.at(i)]++;
    n_neighb_map[verts2.at(i)]++;
  }

  std::vector<EdgeRemovalInfo> edges;
  for (size_t i = 0; i < verts1.size(); ++i) {
    edges.emplace_back(i, n_neighb_map[verts1.at(i)], n_neighb_map[verts2.at(i)], weights.at(i));
  }

  if (Progress::check_abort())
    stop("Interruption");

  if (verbose) {
    std::cout << "Sorting edges..." << std::endl;
  }

  std::sort(edges.begin(), edges.end());

  if (verbose) {
    std::cout << "Filtering edges..." << std::endl;
  }

  std::vector<bool> res(edges.size(), false);
  for (auto const &edge : edges) {
    int i = edge.id;
    if (Progress::check_abort())
      stop("Interruption");

    std::string v1 = verts1.at(i), v2 = verts2.at(i);
    if (v1 == v2) {
      stop("Self-edges are not allowed");
    }

    auto nv1 = n_neighb_map.find(v1);
    auto nv2 = n_neighb_map.find(v2);
    if (nv1->second > min_neighb_per_vertex && nv2->second > min_neighb_per_vertex) {
      int i1 = n_neighb_map.at(v1);
      nv1->second--;
      nv2->second--;
      res.at(edge.id) = true;
    }
  }

  if (verbose) {
    std::cout << "All done!" << std::endl;
  }

  return res;
}