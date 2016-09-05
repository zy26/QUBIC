#include "block.h"
#include "discrete.h"
#include "cluster.h"
#include "version.h"
#include "option.h"
#include "charset.h"
#include "discretize.h"

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstring> // strcmp
class qubic {
  static std::vector<Block> init_qubic(DiscreteArrayListWithSymbols& all, const double c, const double f, std::size_t col_width,
                                       const int o, const Option& option, const CountHelper& count_helper, const bool verbose) {
    if (verbose) fprintf(stdout, "\nQUBIC %s: greedy biclustering\n\n", VER);
    /* ensure enough searching space */
    int SCH_BLOCK = 2 * o;
    /* the file that stores all blocks */
    EdgeList EdgeList(count_helper, verbose);
    /* bi-clustering */
    if (verbose) fprintf(stdout, "Clustering started");
    return cluster(all, EdgeList.get_seeds(), col_width, c, option.cond_, option.area_, option.pvalue_, SCH_BLOCK, o, f, option.filter_1xn_nx1, verbose);
  }

public:
  static std::vector<Block> init_qubic_n(DiscreteArrayListWithSymbols& all, const double c, const double f, std::size_t col_width,
                                         const int o, const Option& option, const bool verbose) {
    return init_qubic(all, c, f, col_width, o, option, CountHelperSaved(all.list, col_width), verbose);
  }

  static std::vector<Block> init_qubic_w(DiscreteArrayListWithSymbols& all, const double c, const double f, std::size_t col_width,
                                         const int o, const Option& option, const bool verbose, const std::vector<std::vector<float>>& weights) {
    return init_qubic(all, c, f, col_width, o, option, WeightedCountHelper(all.list, weights, col_width), verbose);
  }

  static std::vector<Block> init_qubic_e(DiscreteArrayListWithSymbols all, const double c, const std::vector<std::vector<char>>& RowxNumber, const std::vector<std::vector<char>>& NumberxCol) {
    return read_and_solve_blocks(all, c, RowxNumber, NumberxCol);
  }
};

/* Open a file to write or die */
FILE* mustOpenWrite(const char* fileName) {
  if (!strcmp(fileName, "stdin"))  return stdin;
  if (!strcmp(fileName, "stdout")) return stdout;
  FILE* f = fopen(fileName, "w");
  if (f) return f;
  fprintf(stderr, "[Error] Can't open %s to write.", fileName);
  throw - 1;
}

/**************************************************************************/

static void write_imported(const char* stream_nm, const DiscreteArrayList& arr_c, const std::vector<std::string>& genes,
                           const std::vector<std::string>& conds, const DiscreteArray& symbols) {
  FILE* fw = mustOpenWrite(stream_nm);
  fprintf(fw, "o");
  for (std::size_t col = 0; col < conds.size(); col++)
    fprintf(fw, "\t%s", conds[col].c_str());
  fputc('\n', fw);
  for (std::size_t row = 0; row < genes.size(); row++) {
    fprintf(fw, "%s", genes[row].c_str());
    for (std::size_t col = 0; col < conds.size(); col++)
      fprintf(fw, "\t%d", symbols[arr_c[row][col]]);
    fputc('\n', fw);
  }
  fclose(fw);
}

/* Identified clusters are backtraced to the original data, by
* putting the clustered vectors together, identify common column
*/
static void print_bc(FILE* fw, const Block& b, const int& num,
                     const DiscreteArrayList& arr_c, const std::vector<std::string>& genes, const std::vector<std::string>& conds,
                     const Symbols& symbols) {
  int block_rows, block_cols;
  int num_1 = 0, num_2 = 0;
  /* block height (genes) */
  block_rows = b.block_rows();
  /* block_width (conditions) */
  block_cols = b.block_cols();
  fprintf(fw, "BC%03d\tS=%d\tPvalue:%g \n", num, block_rows * block_cols, static_cast<double>(b.pvalue));
  /* fprintf(fw, "BC%03d\tS=%d\tPvalue:%lf \n", num, block_rows * block_cols, (double)b.pvalue); */
  fprintf(fw, " Genes [%d]: ", block_rows);
  for (std::set<int>::iterator it = b.genes_order.begin(); it != b.genes_order.end(); ++it)
    fprintf(fw, "%s ", genes[*it].c_str());
  for (std::set<int>::iterator it = b.genes_reverse.begin(); it != b.genes_reverse.end(); ++it)
    fprintf(fw, "%s ", genes[*it].c_str());
  fprintf(fw, "\n");
  fprintf(fw, " Conds [%d]: ", block_cols);
  for (std::set<int>::iterator it = b.conds.begin(); it != b.conds.end(); ++it)
    fprintf(fw, "%s ", conds[*it].c_str());
  fprintf(fw, "\n");
  /* the complete block data output */
  for (std::set<int>::iterator it = b.genes_order.begin(); it != b.genes_order.end(); ++it) {
    fprintf(fw, "%10s:", genes[*it].c_str());
    for (std::set<int>::iterator jt = b.conds.begin(); jt != b.conds.end(); ++jt) {
      fprintf(fw, "\t%d", symbols[arr_c[*it][*jt]]);
      if (it == b.genes_order.begin()) {
        if (symbols[arr_c[*it][*jt]] == +1) num_1++;
        if (symbols[arr_c[*it][*jt]] == -1) num_2++;
      }
    }
    fputc('\n', fw);
  }
  fputc('\n', fw);
  for (std::set<int>::iterator it = b.genes_reverse.begin(); it != b.genes_reverse.end(); ++it) {
    fprintf(fw, "%10s:", genes[*it].c_str());
    for (std::set<int>::iterator jt = b.conds.begin(); jt != b.conds.end(); ++jt)
      fprintf(fw, "\t%d", symbols[arr_c[*it][*jt]]);
    fputc('\n', fw);
  }
  /*fprintf(stdout, "BC%03d: #of 1 and -1 are:\t%d\t%d\n",num,num_1,num_2);
  fputc('\n', fw);*/
}

void print_params(FILE* fw, const std::size_t col_width, double filter, double TOLERANCE, int RPT_BLOCK) {
  fprintf(fw, "# QUBIC version %.1f output\n", 1.9);
  fprintf(fw, "# \n");
  fprintf(fw, "# Parameters: -k %d -f %.2f -c %.2f -o %d",
          static_cast<unsigned int>(col_width), filter, TOLERANCE, RPT_BLOCK);
  fprintf(fw, "\n\n");
}

/**************************************************************************/

void write_blocks(const std::string& tfile, const std::vector<std::string>& row_names, const std::vector<std::string>& col_names, const double c, const int o, const double filter, int col_width, DiscreteArrayListWithSymbols all, std::vector<Block> output, const bool verbose) {
  FILE* fw = mustOpenWrite((tfile + ".blocks").c_str());
  print_params(fw, col_width, filter, c, o);
  for (std::size_t i = 0; i < output.size(); i++)
    print_bc(fw, output[i], i, all.list, row_names, col_names, all.symbols);
  /* clean up */
  fclose(fw);
  if (verbose) fprintf(stdout, "%d clusters are written to %s\n", static_cast<unsigned int>(output.size()),
                         (tfile + ".blocks").c_str());
}

void write_chars(const std::string& tfile, const std::vector<std::string>& row_names, const std::vector<std::string>& col_names, DiscreteArrayListWithSymbols all, const bool verbose) {
  write_imported((tfile + ".chars").c_str(), all.list, row_names, col_names, all.symbols);
  if (verbose) fprintf(stdout, "Formatted data are written to %s\n", (tfile + ".chars").c_str());
}

std::size_t fix_col_width(const std::vector<std::vector<short>>& x, const int k) {
  return k == 2 ? std::max(x[0].size() / 20, static_cast<std::size_t>(2)) : k;
}

std::vector<Block> r_main(const std::vector<std::vector<short>>& short_matrix, const double c, const int o,
                          const double filter, const int k, const Option& option, const bool verbose) {
  std::size_t col_width = fix_col_width(short_matrix, k);
  if (verbose) fprintf(stdout, "Size of matrix is (%lu, %lu)\n", static_cast<unsigned long>(short_matrix.size()), static_cast<unsigned long>(short_matrix[0].size()));
  DiscreteArrayListWithSymbols all = make_charsets_d(short_matrix, verbose);
  return qubic::init_qubic_n(all, c, filter, col_width, o, option, verbose);
}

std::vector<Block> r_main(const std::vector<std::vector<short>>& short_matrix, const double c, const int o,
                          const double filter, const int k, const Option& option, const bool verbose, const std::vector<std::vector<float>>& weight_matrix) {
  std::size_t col_width = fix_col_width(short_matrix, k);
  if (verbose) fprintf(stdout, "Size of matrix is (%lu, %lu)\n", static_cast<unsigned long>(short_matrix.size()), static_cast<unsigned long>(short_matrix[0].size()));
  DiscreteArrayListWithSymbols all = make_charsets_d(short_matrix, verbose);
  return qubic::init_qubic_w(all, c, filter, col_width, o, option, verbose, weight_matrix);
}

std::vector<Block> r_main(const std::vector<std::vector<short>> &short_matrix, const double c, const bool verbose,
                          const std::vector<std::vector<char>> & RowxNumber, const std::vector<std::vector<char>> & NumberxCol) {
  if (verbose) fprintf(stdout, "Size of matrix is (%lu, %lu)\n", static_cast<unsigned long>(short_matrix.size()), static_cast<unsigned long>(short_matrix[0].size()));
  DiscreteArrayListWithSymbols all = make_charsets_d(short_matrix, verbose);
  return qubic::init_qubic_e(all, c, RowxNumber, NumberxCol);
}

std::vector<Block> main_d(const std::vector<std::vector<short>>& short_matrix, const std::vector<std::string>& row_names,
                          const std::vector<std::string>& col_names, const std::string& tfile,
                          const double c, const int o, const double filter, const int k, const Option& option, const bool verbose) {
  std::size_t col_width = fix_col_width(short_matrix, k);
  DiscreteArrayListWithSymbols all = make_charsets_d(short_matrix, verbose);
  std::vector<Block> output = qubic::init_qubic_n(all, c, filter, col_width, o, option, verbose);
  write_chars(tfile, row_names, col_names, all, verbose);
  write_blocks(tfile, row_names, col_names, c, o, filter, col_width, all, output, verbose);
  return output;
}

void write_rules(const std::string& tfile, const std::vector<std::string>& row_names, std::vector<rule> genes_rules, const bool verbose) {
  FILE* fw = mustOpenWrite((tfile + ".rules").c_str());
  for (std::size_t row = 0; row < row_names.size(); row++) {
    fprintf(fw, "row %s :low=%2.5f, up=%2.5f; %d down-regulated,%d up-regulated\n", row_names[row].c_str(),
            genes_rules[row].lower, genes_rules[row].upper, static_cast<unsigned int>(genes_rules[row].cntl),
            static_cast<unsigned int>(genes_rules[row].cntu));
  }
  fclose(fw);
  if (verbose) fprintf(stdout, "Discretization rules are written to %s\n", (tfile + ".rules").c_str());
}

std::vector<Block> main_c(const std::vector<std::vector<float>>& float_matrix, const std::vector<std::string>& row_names,
                          const std::vector<std::string>& col_names, const std::string& tfile, const short r, const double q,
                          const double c, const int o, const double filter, const int col_width, const Option& option, const bool verbose) {
  std::vector<rule> genes_rules;
  DiscreteArrayList short_matrix = discretize(float_matrix, r, q, genes_rules);
  write_rules(tfile, row_names, genes_rules, verbose);
  return main_d(short_matrix, row_names, col_names, tfile, c, o, filter, col_width, option, verbose);
}
