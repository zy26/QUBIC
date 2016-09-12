#include "qubic.h"
#include "config.h"
#include "discretize.h"
#include "edge_list.h"
#include "fopen_matrix.h"
#include "matrix_float.h"
#include "option.h"
#include <csignal>
#include <algorithm>
#include <iostream>
#include <list>
#include <string>
#include <vector>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

template<typename T>
NumericMatrix from_vector(const std::vector<std::vector<T>> &result) {
  std::size_t nr = result.size();
  std::size_t nc = result[0].size();
  NumericMatrix m(nr, nc);
  for (std::size_t i = 0; i < nr; i++) {
    const std::vector<T> &result_i = result[i];
    if (result_i.size() != nc) stop("QUBIC: incompatible size %d != %d", result_i.size(), nc);
    for (std::size_t j = 0; j < nc; j++) m(i, j) = result_i[j];
  }
  return m;
}

template<typename T, typename TMatrix>
std::vector<std::vector<T>> to_vector(const TMatrix &matrix) {
  auto nc = matrix.ncol();
  auto nr = matrix.nrow();
  std::vector<std::vector<T>> result(nr);
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) result[i].push_back(matrix(i, j));
  }
  return result;
}

template<typename T>
std::vector<std::vector<T>> to_vector2(const arma::sp_mat &matrix) {
  auto nc = matrix.n_cols;
  auto nr = matrix.n_rows;
  std::vector<std::vector<T>> result(nr);
  for (std::size_t i = 0; i < nr; i++) {
    for (std::size_t j = 0; j < nc; j++) result[i].push_back(matrix(i, j));
  }
  return result;
}

NumericMatrix nothing(NumericMatrix matrix) {
  return from_vector<float>(to_vector<float, NumericMatrix>(matrix));
}

List get_list() {
  return List::create();
}

extern "C" void my_function_to_handle_aborts(int signal_number) {
  /*Your code goes here. You can output debugging info.
  If you return from this function, and it was called
  because abort() was called, your program will exit or crash anyway
  (with a dialog box on Windows).
  */
  stop("abort()");
}

List from_blocks(const std::vector<Block> &blocks, const size_t nr, const size_t nc) {
  int number = blocks.size();
  auto x = LogicalMatrix(nr, number);
  auto y = LogicalMatrix(number, nc);
  for (int i = 0; i < number; i++) {
    for (auto it = blocks[i].genes_order.begin(); it != blocks[i].genes_order.end(); ++it) x(*it, i) = true;
    for (auto it = blocks[i].genes_reverse.begin(); it != blocks[i].genes_reverse.end(); ++it) x(*it, i) = true;
    for (auto it = blocks[i].conds.begin(); it != blocks[i].conds.end(); ++it)
      y(i, *it) = true;
  }
  return List::create(
           Named("RowxNumber") = x,
           Named("NumberxCol") = y,
           Named("Number") = blocks.size(),
           Named("info") = get_list());
}

//' @backref src/rcpp_qubic.cpp
// [[Rcpp::export(.qubic_d)]]
List qubic_d(const IntegerMatrix matrix,
             const double c, const int o, const double f, const int k,
             const bool P, const bool S, const bool C,
             const bool verbose) {
  // may treat abort() more friendly, see http://stackoverflow.com/a/3911102
  signal(SIGABRT, &my_function_to_handle_aborts);
  DiscreteArrayList arr_d = to_vector<short, IntegerMatrix>(matrix);
  try {
    std::vector<Block> result = r_main(arr_d, c, o, f, k, Option(P, S, C, true), verbose);
    return from_blocks(result, matrix.nrow(), matrix.ncol());
  }
  catch (double) {
    stop("Something wrong near r_main_d function, maybe out of memory");
  }
  return List::create(); // avoid warning
}

//' @backref src/rcpp_qubic.cpp
// [[Rcpp::export(.qubic_de)]]
List qubic_de(const IntegerMatrix matrix, const double c, 
              const bool verbose, const LogicalMatrix RowxNumber, const LogicalMatrix NumberxCol) {
  // may treat abort() more friendly, see http://stackoverflow.com/a/3911102
  signal(SIGABRT, &my_function_to_handle_aborts);
  try {
    std::vector<Block> result = r_main(to_vector<short, IntegerMatrix>(matrix), c, verbose, to_vector<char, LogicalMatrix>(RowxNumber), to_vector<char, LogicalMatrix>(NumberxCol));
    return from_blocks(result, matrix.nrow(), matrix.ncol());
  }
  catch (double) {
    stop("Something wrong near r_main_d function, maybe out of memory");
  }
  return List::create(); // avoid warning
}

//' @backref src/rcpp_qubic.cpp
// [[Rcpp::export(.qubic_dw)]]
List qubic_dw(const IntegerMatrix matrix,
              const double c, const int o, const double f, const int k,
              const bool P, const bool S, const bool C,
              const bool verbose, const arma::sp_mat weight) {
  // may treat abort() more friendly, see http://stackoverflow.com/a/3911102
  signal(SIGABRT, &my_function_to_handle_aborts);
  try {
    std::vector<Block> result = r_main(to_vector<short, IntegerMatrix>(matrix), c, o, f, k, Option(P, S, C, true), verbose, to_vector2<float>(weight));
    return from_blocks(result, matrix.nrow(), matrix.ncol());
  }
  catch (double) {
    stop("Something wrong near r_main_d function, maybe out of memory");
  }
  return List::create(); // avoid warning
}

//' Create a qualitative discrete matrix for a given gene expression matrix
//'
//' \code{qudiscretize} delivers a discrete matrix. It is useful if we just want to get a discretized matrix.
//'
//' @details
//' \code{qudiscretize} convert a given gene expression matrix to a discrete matrix.
//' It's implimented in C++, providing a increase in speed over the C equivalent.
//'
//' @usage qudiscretize(x, r = 1L, q = 0.06)
//' @inheritParams QUBIC
//'
//' @return A qualitative discrete matrix
//'
//' @name qudiscretize
//'
//' @aliases qudiscretize qdiscretize
//'
//' @examples
//' # Qualitative discretize yeast microarray data
//' data(BicatYeast)
//' qudiscretize(BicatYeast[1:7, 1:5])
//'
//' @seealso \code{\link{QUBIC}} \code{\link{discretize}}
//' @backref src/rcpp_qubic.cpp
// [[Rcpp::export]]
NumericMatrix qudiscretize(const NumericMatrix x, const short r = 1, const double q = 0.06) {
  std::vector<rule> genes_rules;
  auto x1 = to_vector<float, NumericMatrix>(x);
  DiscreteArrayList arr_d = discretize(x1, r, q);
  NumericMatrix result = from_vector(arr_d);
  result.attr("dimnames") = x.attr("dimnames");
  return result;
}
