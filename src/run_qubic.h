#ifndef MAIN_H
#define MAIN_H

#include <string>
#include <vector>
#include "../src/block.h"
#include "../src/matrix_float.h"
#include "../src/option.h"

std::vector<Block> run_qubic_c(const MatrixFloat &matrix, const std::string &tfile = "rQUBIC", const short &r = 1,
  const double &q = 0.06, const double &c = 0.95, const int &o = 100, const double &f = 1, const int &k = 2,
  const Option &option = Option(), const bool verbose = true);
std::vector<Block> run_qubic_d(const Matrix<short> &matrix, const std::string &tfile = "rQUBIC", const double &c = 0.95,
  const int &o = 100, const double &f = 1, const int &k = 2, const Option &option = Option(), const bool verbose = true);


#endif
