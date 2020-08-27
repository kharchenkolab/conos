#include <Rcpp.h>
#include <progress.hpp>

#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix getSumWeightMatrix(const std::vector<double> &weights, const std::vector<int> &row_inds,
                                       const std::vector<int> &col_inds, const std::vector<int> &factor_levels,
                                       bool normalize=true) {
  Rcpp::NumericMatrix m(factor_levels.size(), *std::max_element(factor_levels.begin(), factor_levels.end()));
  for (int i = 0; i < weights.size(); ++i) {
    int cur_row = row_inds.at(i), cur_col = col_inds.at(i);
    int col_fac = factor_levels.at(cur_col) - 1, row_fac = factor_levels.at(cur_row) - 1;
    double edge_weight = weights.at(i);

    m(cur_row, col_fac) += edge_weight;
    m(cur_col, row_fac) += edge_weight;
  }

  if (!normalize)
    return m;

  auto rs = as<Rcpp::NumericVector>(Rcpp::rowSums(m));
  for (int i = 0; i < rs.size(); ++i) {
    m(i, _) = m(i, _) / std::max(rs(i), 1e-10);
  }

  return m;
}

// [[Rcpp::export]]
std::vector<double> adjustWeightsByCellBalancingC(std::vector<double> weights, const std::vector<int> &row_inds, const std::vector<int> &col_inds,
                                                  const std::vector<int> &factor_levels, Rcpp::NumericMatrix dividers) {
  for (int i = 0; i < weights.size(); ++i) {
    int cur_row = row_inds.at(i), cur_col = col_inds.at(i);
    int col_fac = factor_levels.at(cur_col) - 1, row_fac = factor_levels.at(cur_row) - 1;
    weights[i] /= std::sqrt(dividers(cur_row, col_fac) * dividers(cur_col, row_fac));
  }

  return weights;
}
