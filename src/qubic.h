#ifndef QUBIC_H
#define QUBIC_H

#include <string>
#include <vector>

#include "block.h"
#include "option.h"

std::vector<Block> main_c(const std::vector<std::vector<float>> &x, const std::vector<std::string> &row_names,
                          const std::vector<std::string> &col_names, const std::string &tfile, const short r,
                          const double q, const double c, const int o, const double f, const int k, const Option &option, const bool verbose);

std::vector<Block> main_d(const std::vector<std::vector<short>> &x, const std::vector<std::string> &row_names,
                          const std::vector<std::string> &col_names, const std::string &tfile, const double c, const int o, const double f,
                          const int k, const Option &option, const bool verbose);

std::vector<Block> r_main(const std::vector<std::vector<short>> &x,
                          const double c, const int o, const double f, const int k, const Option &option, const bool verbose);

std::vector<Block> r_main(const std::vector<std::vector<short>> &x,
                          const double c, const int o, const double f, const int k, const Option &option, const bool verbose,
                          const std::vector<std::vector<float>> &w);

std::vector<Block> r_main(const std::vector<std::vector<short>> &x,
                          const double c, const bool verbose,
                          const std::vector<std::vector<char>> & RowxNumber, const std::vector<std::vector<char>> & NumberxCol);

#endif
