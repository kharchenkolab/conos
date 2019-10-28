#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <numeric>
#include <Rcpp.h>

#include <progress.hpp>

void trace_time(const std::string &message, bool print_date=false)
{
  std::string format = print_date ? "%m/%d/%Y %H:%M:%S" : "%H:%M:%S";
  time_t ctt = time(nullptr);
  char time_str[100];
  strftime(time_str, sizeof(time_str), format.c_str(), localtime(&ctt));
  Rcpp::Rcout << message << ": " << time_str << "." << std::endl;
}

std::vector<size_t> sortperm(const std::vector<double> &values, bool reverse=false) {
  std::vector<size_t> sorted_ids(values.size());
  std::iota(sorted_ids.begin(), sorted_ids.end(), 0);
  if (reverse) {
    std::sort(sorted_ids.begin(), sorted_ids.end(), [&values](int a, int b){ return values.at(a) > values.at(b); });
  } else {
    std::sort(sorted_ids.begin(), sorted_ids.end(), [&values](int a, int b){ return values.at(a) < values.at(b); });
  }

  return sorted_ids;
}

template<typename T>
std::vector<T> reorder(const std::vector<T> &vec, const std::vector<size_t> indexes, size_t max_size=0)
{
  if (max_size == 0 || max_size > vec.size()) {
    max_size = vec.size();
  }

  auto vec_sorted = std::vector<T>(max_size);
  for (int i = 0; i < max_size; ++i) {
    vec_sorted.at(i) = vec.at(indexes.at(i));
  }

  return vec_sorted;
}

// [[Rcpp::export]]
Rcpp::List as_factor(const std::vector<std::string> &vals) {
  std::unordered_map<std::string, int> levels;
  std::vector<int> int_vals;
  std::vector<std::string> level_vals;

  for (auto const &val : vals) {
    auto const iter = levels.emplace(val, levels.size());
    int_vals.emplace_back(iter.first->second + 1);
    if (iter.second) {
      level_vals.emplace_back(iter.first->first);
    }
  }

  return Rcpp::List::create(Rcpp::_["values"]=int_vals, Rcpp::_["levels"]=level_vals);
}

std::unordered_map<int, double> get_hitting_time_map(const std::vector<int> &adjacent_ids,
                                                     const std::vector<double> &hitting_times) {
  std::unordered_map<int, double> cur_times;
  for (int v2 = 0; v2 < adjacent_ids.size(); ++v2) {
    cur_times[adjacent_ids.at(v2)] = hitting_times.at(v2);
  }

  return cur_times;
}

// Hitting time

class Paths
{
private:
  double _sum_prob;
  double _sum_time;

public:
  Paths(double sum_prob=0, double sum_time=0)
    : _sum_prob(sum_prob)
    , _sum_time(sum_time)
  {}

  void add_path(double probability, int length)
  {
    this->_sum_prob += probability;
    this->_sum_time += probability * length;
  }

  double hitting_distance(double min_prob_threshold=1e-10) const
  {
    double sum_prob = std::max(this->_sum_prob, min_prob_threshold);

    if (sum_prob > 1)
    {
      std::cout << "WARNING: sum of path probabilities is more than 1" << std::endl;
      sum_prob = 1;
    }

    return this->_sum_time / (sum_prob * sum_prob);
  }
};

using path_map_t = std::unordered_map<int, Paths>;

void dfs_hitting_time(const std::vector<std::vector<int>> &adjacency_list,
                      const std::vector<std::vector<double>> &transition_probabilities,
                      int vertex_id, std::vector<bool> &visited, path_map_t &result_paths,
                      double current_prob, int current_length, double min_prob,
                      int min_visited_verts, double min_prob_lower)
{
  if (current_prob < 0 || current_prob > 1 || std::isnan(current_prob)) {
    Rcpp::stop("Wrong current_prob: " + std::to_string(current_prob));
  }

  bool was_visited_before = visited.at(vertex_id);
  if (!was_visited_before)
  {
    auto iter = result_paths.emplace(vertex_id, Paths()).first;
    iter->second.add_path(current_prob, current_length);
  }

  if ((current_prob < min_prob && result_paths.size() >= min_visited_verts) || current_prob < min_prob_lower)
    return;

  visited.at(vertex_id) = true;

  for (int i = 0; i < adjacency_list.at(vertex_id).size(); ++i)
  {
    int neighbor_id = adjacency_list.at(vertex_id).at(i);
    double transition_prob = transition_probabilities.at(vertex_id).at(i);
    if (transition_prob < 0 || transition_prob > 1 || std::isnan(transition_prob)) {
      Rcpp::stop("Wrong transition probability: " + std::to_string(transition_prob));
    }

    dfs_hitting_time(adjacency_list, transition_probabilities, neighbor_id, visited, result_paths,
                     current_prob * transition_prob, current_length + 1, min_prob,
                     min_visited_verts, min_prob_lower);
  }

  visited.at(vertex_id) = was_visited_before;
}

std::pair<std::vector<double>, std::vector<int>> hitting_time_per_neighbor(const std::vector<std::vector<int>> &adjacency_list,
                                                                           const std::vector<std::vector<double>> &transition_probabilities,
                                                                           int start_vertex, double min_prob, int min_visited_verts, double min_prob_lower,
                                                                           int max_adj_num)
{
  if (start_vertex < 0 || start_vertex >= adjacency_list.size())
    Rcpp::stop("Wrong start_vertex index: " + std::to_string(start_vertex));

  path_map_t paths;
  std::vector<bool> visited(adjacency_list.size(), false);
  visited[start_vertex] = true;

  dfs_hitting_time(adjacency_list, transition_probabilities, start_vertex, visited, paths,
                   1.0, 0, min_prob, min_visited_verts, min_prob_lower);

  std::vector<int> indexes;
  std::vector<double> values;
  for (auto const &p : paths)
  {
    indexes.emplace_back(p.first);
    values.emplace_back(p.second.hitting_distance());
  }

  auto sorted_ids = sortperm(values);
  auto values_sorted = reorder(values, sorted_ids, max_adj_num);
  auto indexes_sorted = reorder(indexes, sorted_ids, max_adj_num);

  return std::make_pair(values_sorted, indexes_sorted);
}

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>>
  hitting_time_per_neighbors(const std::vector<std::vector<int>> &adjacency_list,
                             const std::vector<std::vector<double>> &transition_probabilities,
                             int n_verts=0, int n_cores=1, double min_prob = 1e-3,
                             int min_visited_verts=1000, double min_prob_lower=1e-5,
                             int max_adj_num=0, bool verbose=true)
{
#ifdef _OPENMP
  omp_set_num_threads(n_cores);
#endif

  if (n_verts <= 0 || n_verts >= adjacency_list.size()) {
    n_verts = adjacency_list.size();
  }
  // std::vector<Rcpp::NumericVector> res(n_verts);
  std::vector<std::vector<int>> res_idx(n_verts);
  std::vector<std::vector<double>> res_times(n_verts);

  Progress p(n_verts, verbose);

#pragma omp parallel for schedule(dynamic)
  for (int v1 = 0; v1 < n_verts; ++v1) {
    if (Progress::check_abort())
      continue;

    auto const cur_res = hitting_time_per_neighbor(adjacency_list, transition_probabilities, v1,
                                                   min_prob, min_visited_verts, min_prob_lower, max_adj_num);
    p.increment();

#pragma omp critical
{
  res_times.at(v1) = cur_res.first;
  res_idx.at(v1) = cur_res.second;
}
  }

  if (Progress::check_abort())
    Rcpp::stop("Aborted");

  return std::make_pair(res_idx, res_times);
}

// Commute time

Rcpp::List commute_time_per_node(const std::vector<std::vector<int>> &adjacency_list,
                                 const std::vector<std::vector<double>> &hitting_times,
                                 int max_adj_num=0, int n_cores=1, bool verbose=true) {
  if (adjacency_list.size() != hitting_times.size())
    Rcpp::stop("Vectors must have the same length");

#ifdef _OPENMP
  omp_set_num_threads(n_cores);
#endif

  std::vector<std::unordered_map<int, double>> hitting_times_map(adjacency_list.size());

  if (verbose) {
    trace_time("Hashing adjacency list");
  }

  {
    Progress p_hash(adjacency_list.size(), verbose);

#pragma omp parallel for schedule(static)
    for (int v1 = 0; v1 < adjacency_list.size(); ++v1) {
      if (Progress::check_abort())
        continue;

      p_hash.increment();
      hitting_times_map[v1] = get_hitting_time_map(adjacency_list.at(v1), hitting_times.at(v1));
    }

    if (verbose) {
      Rcpp::Rcout << "Done." << std::endl;
    }

    if (Progress::check_abort())
      Rcpp::stop("Aborted");
  }

  std::vector<std::vector<double>> commute_times(adjacency_list.size());
  std::vector<std::vector<int>> commute_time_idx(adjacency_list.size());

  {
    Progress p_dist(hitting_times_map.size(), verbose);

    if (verbose) {
      trace_time("Estimating distances");
    }

#pragma omp parallel for schedule(static)
    for (int v1 = 0; v1 < hitting_times_map.size(); ++v1) {
      if (Progress::check_abort())
        continue;

      p_dist.increment();

      std::vector<double> cur_times;
      std::vector<int> cur_ids;
      for (auto const &p2 : hitting_times_map.at(v1)) {
        auto const &h2_map = hitting_times_map.at(p2.first);
        auto const v1_iter = h2_map.find(v1);
        if (v1_iter == h2_map.end())
          continue;

        cur_times.emplace_back(v1_iter->second + p2.second);
        cur_ids.emplace_back(p2.first);
      }

      auto sorted_ids = sortperm(cur_times);

      commute_times[v1] = reorder(cur_times, sorted_ids, max_adj_num);
      commute_time_idx[v1] = reorder(cur_ids, sorted_ids, max_adj_num);
    }

    if (Progress::check_abort())
      Rcpp::stop("Aborted");
  }

  if (verbose) {
    Rcpp::Rcout << "Done" << std::endl;
  }

#pragma omp barrier

  return Rcpp::List::create(Rcpp::_["idx"]=Rcpp::wrap(commute_time_idx),
                            Rcpp::_["dist"]=Rcpp::wrap(commute_times));
}

// [[Rcpp::export]]
Rcpp::List get_nearest_neighbors(const std::vector<std::vector<int>> &adjacency_list,
                                 const std::vector<std::vector<double>> &transition_probabilities,
                                 int n_verts=0, int n_cores=1, double min_prob = 1e-3,
                                 int min_visited_verts=1000, double min_prob_lower=1e-5,
                                 int max_hitting_nn_num=0, int max_commute_nn_num=0, bool verbose=true) {
  if (verbose) {
    trace_time("Estimating hitting distances");
  }

  auto ht_res = hitting_time_per_neighbors(adjacency_list, transition_probabilities, n_verts,
                                           n_cores, min_prob, min_visited_verts, min_prob_lower,
                                           max_hitting_nn_num, verbose);

  if (verbose) {
    Rcpp::Rcout << "Done." << std::endl;
    trace_time("Estimating commute distances");
  }

  auto res = commute_time_per_node(ht_res.first, ht_res.second, max_commute_nn_num, n_cores, verbose);
  if (verbose) {
    Rcpp::Rcout << "Done." << std::endl;
    trace_time("All done!");
  }

  return res;
}
