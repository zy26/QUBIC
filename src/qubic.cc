#include "qubic.h"
#include "charset.h"
#include "edge_list.h"
#include "version.h"
#include "discretize.h"
#include <cmath> // floor
#include <cstring> // memset
#include <cstddef> // size_t
#include <algorithm>
#include <list>

namespace internal
{
  static void seed_update(const DiscreteArray& s, std::vector<std::vector<bits16>>& profile)
  {
    for (size_t i = 0; i < s.size(); i++)
      profile[i][s[i]]++;
  }

  /* scan through all columns and identify the set within threshold,
     * "fuzziness" of the block is controlled by TOLERANCE (-c)
     */
  static void scan_block(const DiscreteArrayList& arr_c, const Symbols& symbols, const std::vector<int>& gene_set,
                         Block& b, std::vector<std::vector<bits16>>& profile, double TOLERANCE)
  {
    size_t i, j;
    size_t block_rows, cur_rows;
    block_rows = cur_rows = gene_set.size();
    size_t k;
    for (j = 0; j < profile.size(); j++)
      for (k = 0; k < profile[j].size(); k++)
        profile[j][k] = 0;
    for (j = 0; j < cur_rows; j++)
      seed_update(arr_c[gene_set[j]], profile);
    int btolerance = static_cast<int>(ceil(TOLERANCE * block_rows));
    for (j = 0; j < profile.size(); j++)
    {
      /* See if this column satisfies tolerance */
      /* here i start from 1 because symbols[0]=0 */
      for (i = 1; i < symbols.size(); i++)
      {
        if ((profile[j][i] >= btolerance))
        {
          b.conds.insert(j);
          break;
        }
      }
    }
  }

  /*************************************************************************/

  static void update_colcand(std::list<size_t>& colcand, const DiscreteArray& g1, const DiscreteArray& g2)
  {
    std::list<size_t>::iterator it = colcand.begin();
    while (it != colcand.end())
    {
      if ((g1[*it] != g2[*it]))
        colcand.erase(it++); // alternatively, it = colcand.erase(it);
      else
        ++it;
    }
  }

  /*calculate the weight of the edge with two vertices g1 and g2*/
  static int intersect_row(const std::list<size_t>& colcand, const DiscreteArray& g1, const DiscreteArray& g2)
  {
    int cnt = 0;
    for (auto it = colcand.begin(); it != colcand.end(); ++it)
      if ((g1[*it] == g2[*it]) && (g1[*it] != 0))
        cnt++;
    return cnt;
  }

  /*calculate the negative correlation between g1 and g2*/
  inline int reverse_row(const std::list<size_t>& colcand, const DiscreteArray& g1, const DiscreteArray& g2,
                         const std::vector<discrete>& symbols)
  {
    int cnt = 0;
    for (auto it = colcand.begin(); it != colcand.end(); ++it) if ((symbols[g1[*it]] == -symbols[g2[*it]])) cnt++;
    return cnt;
  }

  /* calculate the coverage of any row to the current consensus
    * cnt = # of valid consensus columns
    */
  static int seed_current_modify(const DiscreteArray& s, std::list<size_t>& colcand, const int components,
                                 std::vector<std::vector<bits16>>& profile, double TOLERANCE)
  {
    size_t n;
    int threshold = static_cast<int>(ceil(components * TOLERANCE));
    discrete ss;
    int cnt = 0;
    for (size_t i = 0; i < profile.size(); i++)
    {
      size_t flag = 0;
      ss = s[i];
      for (size_t k = 1; k < profile[i].size(); k++)
      {
        n = profile[i][k];
        if (static_cast<discrete>(k) == ss) n++;
        if (static_cast<int>(n) >= threshold)
        {
          flag = k;
          break;
        }
      }
      if (flag)
      {
        cnt++;
        colcand.push_back(i);
      }
    }
    return cnt;
  }

  /*check whether current edge can be treat as a seed*/
  static bool check_seed(const Edge* e, const std::vector<Block>& bb, size_t rows)
  {
    int block_id = bb.size();
    int i, b1, b2, b3;
    bool fg = false;
    b1 = b2 = -1;
    for (i = 0; i < block_id; i++)
      if ((bb[i].contains(e->gene_one)) && (bb[i].contains(e->gene_two)))
        return false;
    std::vector<int> profiles(rows, 0);
    fg = false;
    for (i = 0; i < block_id; i++)
      if (bb[i].contains(e->gene_one))
      {
        fg = true;
        break;
      }
    if (fg)
      b1 = i;
    fg = false;
    for (i = 0; i < block_id; i++)
      if (bb[i].contains(e->gene_two))
      {
        fg = true;
        break;
      }
    if (fg)
      b2 = i;
    if ((b1 == -1) || (b2 == -1))
      return true;
    else
    {
      for (std::set<int>::iterator it = bb[b1].genes_order.begin(); it != bb[b1].genes_order.end(); ++it)
        profiles[*it]++;
      for (std::set<int>::iterator it = bb[b1].genes_reverse.begin(); it != bb[b1].genes_reverse.end(); ++it)
        profiles[*it]++;
      for (std::set<int>::iterator it = bb[b2].genes_order.begin(); it != bb[b2].genes_order.end(); ++it)
        profiles[*it]++;
      for (std::set<int>::iterator it = bb[b2].genes_reverse.begin(); it != bb[b2].genes_reverse.end(); ++it)
        profiles[*it]++;
      for (size_t index = 0; index < rows; index++)
        if (profiles[index] > 1) return false;
      b3 = std::max(bb[b1].block_cols(), bb[b2].block_cols());
      return !(e->score < b3/* (bb[b1]->block_cols + bb[b2]->block_cols) / 2*/);
    }
  }

  static long double get_pvalue(const continuous& a, const int& b)
  {
    long double pvalue = 0;
    long double poisson = 1.0 / exp(a);
    for (int i = 0; i < b + 300; i++)
    {
      if (i > (b - 1)) pvalue += poisson;
      else poisson *= a / (i + 1.0);
    }
    return pvalue;
  }

  static void block_init(const DiscreteArrayList& arr_c, Block& b,
                         std::vector<int>& genes, std::vector<int>& scores,
                         std::vector<bool>& candidates, const int& cand_threshold,
                         size_t& components, std::vector<long double>& pvalues, bool IS_cond, size_t COL_WIDTH, bool IS_area)
  {
    size_t rows = arr_c.size();
    size_t cols = arr_c[0].size();
    int score, top;
    int cnt = 0, cnt_all = 0, pid = 0;
    continuous cnt_ave = 0, row_all = static_cast<continuous>(rows);
    long double pvalue;
    int max_cnt, max_i;
    std::vector<int> arr_rows(rows), arr_rows_b(rows);
    std::list<size_t> colcand;
    DiscreteArray g1, g2;
    g1 = arr_c[genes[0]];
    g2 = arr_c[genes[1]];
    for (size_t i = 0; i < cols; i++)
      if ((g1[i] == g2[i]) && (g1[i] != 0))
        colcand.push_back(i);
    for (size_t i = 0; i < rows; i++)
    {
      arr_rows[i] = intersect_row(colcand, arr_c[genes[0]], arr_c[i]);
      arr_rows_b[i] = arr_rows[i];
    }
    /*we just get the largest 100 rows when we initial a bicluster because we believe that
       * the 100 rows can characterize the structure of the bicluster
       * btw, it can reduce the time complexity*/
    if (rows > 100)
    {
      std::sort(arr_rows_b.begin(), arr_rows_b.end());
      top = arr_rows_b[rows - 100];
      for (size_t i = 0; i < rows; i++)
        if (arr_rows[i] < top)
          candidates[i] = false;
    }
    /*calculate the condition low bound for current seed*/
    int cutoff = static_cast<int>(0.05 * rows);
    b.cond_low_bound = arr_rows_b[rows - cutoff - 1];
    while (components < rows)
    {
      max_cnt = -1;
      max_i = -1;
      components++;
      cnt_all = 0;
      cnt_ave = 0;
      /******************************************************/
      /*add a function of controlling the bicluster by pvalue*/
      /******************************************************/
      for (size_t i = 0; i < rows; i++)
      {
        if (!candidates[i]) continue;
        cnt = intersect_row(colcand, arr_c[genes[0]], arr_c[i]);
        cnt_all += cnt;
        if (cnt < cand_threshold)
          candidates[i] = false;
        if (cnt > max_cnt)
        {
          max_cnt = cnt;
          max_i = i;
        }
      }
      cnt_ave = cnt_all / row_all;
      pvalue = get_pvalue(cnt_ave, max_cnt);
      if (IS_cond)
      {
        if (max_cnt < static_cast<int>(COL_WIDTH) || max_i < 0 || max_cnt < b.cond_low_bound) break;
      }
      else
      {
        if (max_cnt < static_cast<int>(COL_WIDTH) || max_i < 0) break;
      }
      if (IS_area) score = components * max_cnt;
      else score = std::min(static_cast<int>(components), max_cnt);
      if (score > b.score) b.score = score;
      if (pvalue < b.pvalue) b.pvalue = pvalue;
      genes.push_back(max_i);
      scores.push_back(score);
      pvalues[pid++] = pvalue;
      update_colcand(colcand, arr_c[genes[0]], arr_c[max_i]);
      candidates[max_i] = false;
    }
    /*be sure to free a pointer when you finish using it*/
  }

  /* compare function for qsort, descending by score */
  static bool block_cmpr(const Block& a, const Block& b)
  {
    return a.score > b.score;
  }

  /************************************************************************/
  static std::vector<Block> report_blocks(std::vector<Block> bb, int RPT_BLOCK, double FILTER)
  {
    std::vector<Block> output;
    int num = bb.size();
    std::sort(bb.begin(), bb.end(), block_cmpr);
    int i, j, k;
    /*MIN MAX et al functions can be accessed in struct.h*/
    int n = std::min(num, RPT_BLOCK);
    bool flag;
    double cur_rows, cur_cols;
    double inter_rows, inter_cols;
    /*double proportion;*/
    /* the major post-processing here, filter overlapping blocks*/
    i = 0;
    j = 0;
    while (i < num && j < n)
    {
      Block& b_ptr = bb[i];
      cur_rows = b_ptr.block_rows();
      cur_cols = b_ptr.block_cols();
      flag = true;
      k = 0;
      while (k < j)
      {
        inter_rows = count_intersect(output[k].genes_order, b_ptr.genes_order) +
          count_intersect(output[k].genes_order, b_ptr.genes_reverse) +
          count_intersect(output[k].genes_reverse, b_ptr.genes_order) +
          count_intersect(output[k].genes_reverse, b_ptr.genes_reverse);
        inter_cols = count_intersect(output[k].conds, b_ptr.conds);
        if (inter_rows * inter_cols > FILTER * cur_rows * cur_cols)
        {
          flag = false;
          break;
        }
        k++;
      }
      i++;
      if (flag)
      {
        output.push_back(b_ptr);
        j++;
      }
    }
    return output;
  }

  /************************************************************************/
  static std::vector<Block> cluster(const DiscreteArrayListWithSymbols& all, const std::vector<Edge *>& el,
                                    size_t COL_WIDTH, double TOLERANCE, bool IS_cond, bool IS_area,
                                    bool IS_pvalue, size_t SCH_BLOCK, int RPT_BLOCK, double FILTER, double f, bool verbose)
  {
    std::vector<Block> bb;
    size_t rows = all.list.size();
    size_t cols = all.list[0].size();
    std::vector<long double> pvalues(rows);
    std::vector<bool> candidates(rows);
    std::set<int> allincluster;
    for (std::vector<Edge *>::const_iterator it = el.begin(); it != el.end(); ++it)
    {
      const Edge* e = *it;
      /* check if both genes already enumerated in previous blocks */
      bool flag = true;
      /* speed up the program if the rows bigger than 200 */
      if (rows > 250)
      {
        if (allincluster.find(e->gene_one) != allincluster.end()
          && allincluster.find(e->gene_two) != allincluster.end())
          flag = false;
      }
      else flag = check_seed(e, bb, rows);
      if (!flag) continue;
      std::vector<std::vector<bits16>> profile(cols, std::vector<bits16>(all.symbols.size(), 0));
      /*you must allocate a struct if you want to use the pointers related to it*/
      Block b;
      /*initial the b->score*/
      b.score = std::min(2, e->score);
      /*initial the b->pvalue*/
      b.pvalue = 1;
      /* initialize the stacks genes and scores */
      std::vector<int> genes_order, genes_reverse, scores;
      genes_order.reserve(rows);
      genes_reverse.reserve(rows);
      scores.reserve(rows);
      genes_order.push_back(e->gene_one);
      genes_order.push_back(e->gene_two);
      scores.push_back(1);
      scores.push_back(b.score);
      /* branch-and-cut condition for seed expansion */
      int cand_threshold = static_cast<int>(floor(COL_WIDTH * TOLERANCE));
      if (cand_threshold < 2)
        cand_threshold = 2;
      /* maintain a candidate list to avoid looping through all rows */
      for (size_t j = 0; j < rows; j++)
        candidates[j] = true;
      candidates[e->gene_one] = candidates[e->gene_two] = false;
      size_t components = 2;
      /* expansion step, generate a bicluster without noise */
      block_init(all.list, b, genes_order, scores, candidates, cand_threshold, components, pvalues, IS_cond, COL_WIDTH,
                 IS_area);
      /* track back to find the genes by which we get the best score*/
      size_t k;
      for (k = 0; k < components; k++)
      {
        if (IS_pvalue)
          if ((pvalues[k] == b.pvalue) && (k >= 2) && (scores[k] != scores[k + 1])) break;
        if ((scores[k] == b.score) && ((k + 1 == scores.size()) || (scores[k + 1] != b.score))) break;
      }
      components = k + 1;
      std::fill(candidates.begin(), candidates.end(), true);
      for (size_t ki = 0; ki < components - 1; ki++)
      {
        seed_update(all.list[genes_order[ki]], profile);
        candidates[genes_order[ki]] = false;
      }
      candidates[genes_order[k]] = false;
      genes_order.resize(k + 1);
      std::list<size_t> colcand;
      /* add columns satisfy the conservative r */
      int cnt = seed_current_modify(all.list[genes_order[k]], colcand, components, profile, TOLERANCE);
      /* add some new possible genes */
      int m_cnt;
      for (size_t ki = 0; ki < rows; ki++)
      {
        m_cnt = intersect_row(colcand, all.list[genes_order[0]], all.list[ki]);
        if (candidates[ki] && (m_cnt >= floor(cnt * TOLERANCE)))
        {
          genes_order.push_back(ki);
          components++;
          candidates[ki] = false;
        }
      }
      /* add genes that negative regulated to the consensus */
      for (size_t ki = 0; ki < rows; ki++)
      {
        m_cnt = reverse_row(colcand, all.list[genes_order[0]], all.list[ki], all.symbols);
        if (candidates[ki] && (m_cnt >= floor(cnt * TOLERANCE)))
        {
          genes_reverse.push_back(ki);
          components++;
          candidates[ki] = false;
        }
      }
      /* store gene arrays inside block */
      scan_block(all.list, all.symbols, genes_order, b, profile, TOLERANCE);
      if (b.block_cols() == 0) continue;
      //b.block_rows = components;
      if (IS_pvalue) b.score = static_cast<int>(-(100 * log(b.pvalue)));
      else b.score = components * b.block_cols();
      for (std::vector<int>::iterator iterator = genes_order.begin(); iterator != genes_order.end(); ++iterator)
      {
        b.genes_order.insert(*iterator);
        allincluster.insert(*iterator);
      }
      for (std::vector<int>::iterator iterator = genes_reverse.begin(); iterator != genes_reverse.end(); ++iterator)
      {
        b.genes_reverse.insert(*iterator);
        allincluster.insert(*iterator);
      }

      if (f && (b.block_cols() <= 1 || b.block_rows() <= 1)) continue;

      /*save the current block b to the block list bb so that we can sort the blocks by their score*/
      bb.push_back(b);
      /* reaching the results number limit */
      if (bb.size() == SCH_BLOCK) break;
      if (verbose) fputc('.', stdout);
    }
    if (verbose) fprintf(stdout, "\n");
    return report_blocks(bb, RPT_BLOCK, FILTER);
  }
}

class qubic
{
public:
  static std::vector<Block> init_qubic(DiscreteArrayListWithSymbols& all, const double c, const double f, size_t k,
                                       const int o, const Option& option, const bool verbose)
  {
    if (verbose) fprintf(stdout, "\nQUBIC %s: greedy biclustering\n\n", VER);
    /* ensure enough searching space */
    int SCH_BLOCK = 2 * o;
    /* case 's': po.IS_SWITCH = true; */
    bool IS_pvalue = option.p;
    bool IS_area = option.s;
    bool IS_cond = option.c;
    /* the file that stores all blocks */
    EdgeList EdgeList(all.list, k, verbose);
    /* bi-clustering */
    if (verbose) fprintf(stdout, "Clustering started");
    return internal::cluster(all, EdgeList.get_edge_list(), k, c, IS_cond, IS_area, IS_pvalue, SCH_BLOCK, o, f, option.filter_1xn_nx1, verbose);
  }
};

/**************************************************************************/
/* file-related operations */

/* Strings */
/* strcmp: a zero value indicates that both strings are equal.
* a value greater than zero indicates that the first character that does not match has a greater value in str1 than in str2;
* And a value less than zero indicates the opposite.
*/
#define sameString(a, b) (strcmp((a), (b))==0)

/* Returns TRUE if two strings are same */

/* Open a file to write or die */
FILE* mustOpenWrite(const char* fileName)
{
  FILE* f;
  if (sameString(fileName, "stdin")) return stdin;
  if (sameString(fileName, "stdout")) return stdout;
  if ((f = fopen(fileName, "w")) == NULL)
  {
    fprintf(stderr, "[Error] Can't open %s to write.", fileName);
    throw - 1;
  }
  return f;
}

/**************************************************************************/

static void write_imported(const char* stream_nm, const DiscreteArrayList& arr_c, const std::vector<std::string>& genes,
                           const std::vector<std::string>& conds, const std::vector<discrete>& symbols)
{
  FILE* fw = mustOpenWrite(stream_nm);
  fprintf(fw, "o");
  for (size_t col = 0; col < conds.size(); col++)
    fprintf(fw, "\t%s", conds[col].c_str());
  fputc('\n', fw);
  for (size_t row = 0; row < genes.size(); row++)
  {
    fprintf(fw, "%s", genes[row].c_str());
    for (size_t col = 0; col < conds.size(); col++)
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
                     const std::vector<discrete>& symbols)
{
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
  for (std::set<int>::iterator it = b.genes_order.begin(); it != b.genes_order.end(); ++it)
  {
    fprintf(fw, "%10s:", genes[*it].c_str());
    for (std::set<int>::iterator jt = b.conds.begin(); jt != b.conds.end(); ++jt)
    {
      fprintf(fw, "\t%d", symbols[arr_c[*it][*jt]]);
      if (it == b.genes_order.begin())
      {
        if (symbols[arr_c[*it][*jt]] == 1) num_1++;
        if (symbols[arr_c[*it][*jt]] == -1) num_2++;
      }
    }
    fputc('\n', fw);
  }
  fputc('\n', fw);
  for (std::set<int>::iterator it = b.genes_reverse.begin(); it != b.genes_reverse.end(); ++it)
  {
    fprintf(fw, "%10s:", genes[*it].c_str());
    for (std::set<int>::iterator jt = b.conds.begin(); jt != b.conds.end(); ++jt)
      fprintf(fw, "\t%d", symbols[arr_c[*it][*jt]]);
    fputc('\n', fw);
  }
  /*fprintf(stdout, "BC%03d: #of 1 and -1 are:\t%d\t%d\n",num,num_1,num_2);
  fputc('\n', fw);*/
}

void print_params(FILE* fw, const size_t COL_WIDTH, double FILTER, double TOLERANCE, int RPT_BLOCK)
{
  fprintf(fw, "# QUBIC version %.1f output\n", 1.9);
  fprintf(fw, "# \n");
  fprintf(fw, "# Parameters: -k %d -f %.2f -c %.2f -o %d",
          static_cast<unsigned int>(COL_WIDTH), FILTER, TOLERANCE, RPT_BLOCK);
  fprintf(fw, "\n\n");
}

/**************************************************************************/

std::vector<Block> main_d(const std::vector<std::vector<short>>& x, const std::vector<std::string>& row_names,
                          const std::vector<std::string>& col_names, const std::string& tfile,
                          const double c, const int o, const double f, const int k, const Option& option, const bool verbose)
{
  int k1;
  if (k == 2) k1 = std::max(x[0].size() / 20, static_cast<size_t>(2)); else k1 = k;
  DiscreteArrayListWithSymbols all = make_charsets_d(x, verbose);
  std::vector<Block> output = qubic::init_qubic(all, c, f, k1, o, option, verbose);
  write_imported((tfile + ".chars").c_str(), all.list, row_names, col_names, all.symbols);
  if (verbose) fprintf(stdout, "Formatted data are written to %s\n", (tfile + ".chars").c_str());
  {
    FILE* fw = mustOpenWrite((tfile + ".blocks").c_str());
    print_params(fw, k1, f, c, o);
    for (size_t i = 0; i < output.size(); i++)
      print_bc(fw, output[i], i, all.list, row_names, col_names, all.symbols);
    /* clean up */
    fclose(fw);
    if (verbose)
      fprintf(stdout, "%d clusters are written to %s\n", static_cast<unsigned int>(output.size()),
                    (tfile + ".blocks").c_str());
  }
  return output;
}

std::vector<Block> main_c(const std::vector<std::vector<float>>& x, const std::vector<std::string>& row_names,
                          const std::vector<std::string>& col_names, const std::string& tfile, const short r, const double q,
                          const double c, const int o, const double f, const int k, const Option& option, const bool verbose)
{
  std::vector<rule> genes_rules;
  std::vector<std::vector<discrete>> arr_d = discretize(x, q, r, genes_rules);
  {
    FILE* fw = mustOpenWrite((tfile + ".rules").c_str());
    for (size_t row = 0; row < row_names.size(); row++)
    {
      fprintf(fw, "row %s :low=%2.5f, up=%2.5f; %d down-regulated,%d up-regulated\n", row_names[row].c_str(),
              genes_rules[row].lower, genes_rules[row].upper, static_cast<unsigned int>(genes_rules[row].cntl),
              static_cast<unsigned int>(genes_rules[row].cntu));
    }
    fclose(fw);
    if (verbose) fprintf(stdout, "Discretization rules are written to %s\n", (tfile + ".rules").c_str());
  }
  return main_d(arr_d, row_names, col_names, tfile, c, o, f, k, option, verbose);
}

std::vector<Block> r_main_d(const std::vector<std::vector<short>>& x, const double c, const int o,
                            const double f, const int k, const Option& option, const bool verbose)
{
  if (verbose) fprintf(stdout, "Size of matrix is (%lu, %lu)\n", static_cast<unsigned long>(x.size()), static_cast<unsigned long>(x[0].size()));
  DiscreteArrayListWithSymbols all = make_charsets_d(x, verbose);
  return qubic::init_qubic(all, c, f, k, o, option, verbose);
}

std::vector<Block> r_main_c(const std::vector<std::vector<float>>& x, const short r, const double q,
                            const double c, const int o, const double f, const int k, const Option& option, const bool verbose)
{
  std::vector<rule> genes_rules;
  std::vector<std::vector<discrete>> arr_d = discretize(x, q, r, genes_rules);
  return r_main_d(arr_d, c, o, f, k, option, verbose);
}
