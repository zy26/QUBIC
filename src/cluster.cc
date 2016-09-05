#include "cluster.h"
#include "edge_list.h"
#include "discretize.h"
#include <cmath> // floor, ceil
#include <cstddef> // size_t
#include <algorithm>
#include <list>
#include <vector>

namespace internal {
static void seed_update(const DiscreteArray& s, std::vector<std::vector<bits16>>& profile) {
  for (std::size_t i = 0; i < s.size(); i++)
    profile[i][s[i]]++;
}

/* scan through all columns and identify the set within threshold,
* "fuzziness" of the block is controlled by TOLERANCE (-c)
*/
void scan_block(const DiscreteArrayList& arr_c, const Symbols& symbols, const std::vector<int>& gene_set,
                Block& b, double TOLERANCE) {
  assert(symbols.size() <= std::numeric_limits<unsigned int>::max());
  assert(arr_c[0].size() <= std::numeric_limits<unsigned int>::max());
  unsigned int size = static_cast<unsigned int>(arr_c[0].size());
  unsigned int count = static_cast<unsigned int>(symbols.size());
  std::size_t block_rows;
  block_rows = gene_set.size();
  std::vector<std::vector<bits16>> profile((size), std::vector<bits16>(count, 0));
  for (std::size_t j = 0; j < gene_set.size(); j++)
    seed_update(arr_c[gene_set[j]], profile);
  int btolerance = static_cast<int>(std::ceil(TOLERANCE * block_rows));
  for (unsigned int j = 0; j < size; j++) {
    /* See if this column satisfies tolerance */
    /* here i start from 1 because symbols[0]=0 */
    for (unsigned int i = 1; i < count; i++) {
      if (profile[j][i] >= btolerance) {
        b.conds.insert(j);
        break;
      }
    }
  }
}

/*************************************************************************/

static void update_colcand(std::list<std::size_t>& colcand, const DiscreteArray& g1, const DiscreteArray& g2) {
  std::list<std::size_t>::iterator it = colcand.begin();
  while (it != colcand.end()) {
    if (g1[*it] != g2[*it]) colcand.erase(it++); // alternatively, it = colcand.erase(it);
    else ++it;
  }
}

/*calculate the weight of the edge with two vertices g1 and g2*/
static int intersect_row(const std::list<std::size_t>& colcand, const DiscreteArray& g1, const DiscreteArray& g2) {
  int cnt = 0;
  for (auto it = colcand.begin(); it != colcand.end(); ++it) if (g1[*it] == g2[*it] && g1[*it] != 0) cnt++;
  return cnt;
}

/*calculate the negative correlation between g1 and g2*/
inline int reverse_row(const std::list<std::size_t>& colcand, const DiscreteArray& g1, const DiscreteArray& g2,
                       const Symbols& symbols) {
  int cnt = 0;
  for (auto it = colcand.begin(); it != colcand.end(); ++it) if (symbols[g1[*it]] == -symbols[g2[*it]]) cnt++;
  return cnt;
}

/* calculate the coverage of any row to the current consensus
* cnt = # of valid consensus columns
*/
static int seed_current_modify(const DiscreteArray& s, std::list<std::size_t>& colcand, const int components,
                               std::vector<std::vector<bits16>>& profile, double TOLERANCE) {
  std::size_t n;
  int threshold = static_cast<int>(std::ceil(components * TOLERANCE));
  int cnt = 0;
  for (std::size_t i = 0; i < profile.size(); i++) {
    std::size_t flag = 0;
    short ss = s[i];
    for (std::size_t k = 1; k < profile[i].size(); k++) {
      n = profile[i][k];
      if (static_cast<short>(k) == ss) n++;
      if (static_cast<int>(n) >= threshold) {
        flag = k;
        break;
      }
    }
    if (flag) {
      cnt++;
      colcand.push_back(i);
    }
  }
  return cnt;
}

/*check whether current edge can be treat as a seed*/
static bool check_seed(const Edge* e, const std::vector<Block>& bb, std::size_t rows) {
  std::size_t block_id = bb.size();
  std::size_t b1, b2;
  std::size_t b3;
  for (std::size_t i = 0; i < block_id; i++)
    if (bb[i].contains(e->gene_one) && bb[i].contains(e->gene_two))
      return false;
  std::vector<int> profiles(rows, 0);
  bool fg = false;
  for (std::size_t i = 0; i < block_id; i++)
    if (bb[i].contains(e->gene_one)) {
      fg = true;
      b1 = i;
      break;
    }
  if (!fg) return true;
  fg = false;
  for (std::size_t i = 0; i < block_id; i++)
    if (bb[i].contains(e->gene_two)) {
      fg = true;
      b2 = i;
      break;
    }
  if (!fg) return true;
  for (std::set<int>::iterator it = bb[b1].genes_order.begin(); it != bb[b1].genes_order.end(); ++it)
    profiles[*it]++;
  for (std::set<int>::iterator it = bb[b1].genes_reverse.begin(); it != bb[b1].genes_reverse.end(); ++it)
    profiles[*it]++;
  for (std::set<int>::iterator it = bb[b2].genes_order.begin(); it != bb[b2].genes_order.end(); ++it)
    profiles[*it]++;
  for (std::set<int>::iterator it = bb[b2].genes_reverse.begin(); it != bb[b2].genes_reverse.end(); ++it)
    profiles[*it]++;
  for (std::size_t index = 0; index < rows; index++)
    if (profiles[index] > 1) return false;
  b3 = std::max(bb[b1].block_cols(), bb[b2].block_cols());
  return e->score - b3 >= 0;
}

static long double get_pvalue(const continuous& a, const int& b) {
  long double pvalue = 0;
  long double poisson = 1.0 / exp(a);
  for (int i = 0; i < b + 300; i++) {
    if (i > b - 1) pvalue += poisson;
    else poisson *= a / (i + 1.0);
  }
  return pvalue;
}

static void block_init(const DiscreteArrayList& arr_c, Block& b,
                       std::vector<int>& genes, std::vector<int>& scores,
                       std::vector<char>& candidates, const int& cand_threshold,
                       std::size_t& components, std::vector<long double>& pvalues, bool IS_cond, std::size_t COL_WIDTH, bool IS_area) {
  std::size_t rows = arr_c.size();
  std::size_t cols = arr_c[0].size();
  int score, top;
  continuous cnt_ave = 0, row_all = static_cast<continuous>(rows);
  long double pvalue;
  int max_cnt, max_i;
  std::vector<int> arr_rows(rows), arr_rows_b(rows);
  std::list<std::size_t> colcand;
  DiscreteArray g1, g2;
  g1 = arr_c[genes[0]];
  g2 = arr_c[genes[1]];
  for (std::size_t i = 0; i < cols; i++)
    if (g1[i] == g2[i] && g1[i] != 0)
      colcand.push_back(i);
  for (std::size_t i = 0; i < rows; i++) {
    arr_rows[i] = intersect_row(colcand, arr_c[genes[0]], arr_c[i]);
    arr_rows_b[i] = arr_rows[i];
  }
  /*we just get the largest 100 rows when we initial a bicluster because we believe that
  * the 100 rows can characterize the structure of the bicluster
  * btw, it can reduce the time complexity*/
  if (rows > 100) {
    std::sort(arr_rows_b.begin(), arr_rows_b.end());
    top = arr_rows_b[rows - 100];
    for (std::size_t i = 0; i < rows; i++)
      if (arr_rows[i] < top)
        candidates[i] = false;
  }
  /*calculate the condition low bound for current seed*/
  int cutoff = static_cast<int>(0.05 * rows);
  b.cond_low_bound = arr_rows_b[rows - cutoff - 1];
  while (components < rows) {
    max_cnt = -1;
    max_i = -1;
    components++;
    int cnt_all = 0;
    cnt_ave = 0;
    /******************************************************/
    /*add a function of controlling the bicluster by pvalue*/
    /******************************************************/
    for (std::size_t i = 0; i < rows; i++) {
      if (!candidates[i]) continue;
      int cnt = intersect_row(colcand, arr_c[genes[0]], arr_c[i]); // TODO (Qin): Yu want to find out why not reverse intersect rows are considered.
      cnt_all += cnt;
      if (cnt < cand_threshold)
        candidates[i] = false;
      if (cnt > max_cnt) {
        max_cnt = cnt;
        max_i = i;
      }
    }
    cnt_ave = cnt_all / row_all;
    if (max_i < 0) break;
    if (max_cnt < static_cast<int>(COL_WIDTH)) break;
    if (IS_cond && max_cnt < b.cond_low_bound) break;
    if (IS_area) score = components * max_cnt;
    else score = std::min(static_cast<int>(components), max_cnt);
    if (score > b.score) b.score = score;
    pvalue = get_pvalue(cnt_ave, max_cnt);
    if (pvalue < b.pvalue) b.pvalue = pvalue;
    genes.push_back(max_i);
    scores.push_back(score);
    pvalues.push_back(pvalue);
    update_colcand(colcand, arr_c[genes[0]], arr_c[max_i]);
    candidates[max_i] = false;
  }
  /*be sure to free a pointer when you finish using it*/
}

/* compare function for qsort, descending by score */
static bool block_cmpr(const Block& a, const Block& b) {
  return a.score > b.score;
}

/************************************************************************/
static std::vector<Block> report_blocks(std::vector<Block> bb, int RPT_BLOCK, double FILTER) {
  std::vector<Block> output;
  int num = bb.size();
  std::sort(bb.begin(), bb.end(), block_cmpr);
  /*MIN MAX et al functions can be accessed in struct.h*/
  int n = std::min(num, RPT_BLOCK);
  /*double proportion;*/
  /* the major post-processing here, filter overlapping blocks*/
  int i = 0;
  int j = 0;
  while (i < num && j < n) {
    Block& b_ptr = bb[i];
    double cur_rows = b_ptr.block_rows();
    double cur_cols = b_ptr.block_cols();
    bool flag = true;
    int k = 0;
    while (k < j) {
      double inter_rows = count_intersect(output[k].genes_order, b_ptr.genes_order) +
                          count_intersect(output[k].genes_order, b_ptr.genes_reverse) +
                          count_intersect(output[k].genes_reverse, b_ptr.genes_order) +
                          count_intersect(output[k].genes_reverse, b_ptr.genes_reverse);
      double inter_cols = count_intersect(output[k].conds, b_ptr.conds);
      if (inter_rows * inter_cols > FILTER * cur_rows * cur_cols) {
        flag = false;
        break;
      }
      k++;
    }
    i++;
    if (flag) {
      output.push_back(b_ptr);
      j++;
    }
  }
  return output;
}
}

int add_intersect(const DiscreteArrayListWithSymbols& all, std::vector<int>& genes_order, std::vector<char>& candidates,
                  const std::list<std::size_t>& colcand, const DiscreteArray& first, double threadshold) {
  std::size_t rows = all.list.size();
  int count = 0;
  for (std::size_t i = 0; i < rows; i++) {
    if (!candidates[i]) continue;
    int m_cnt = internal::intersect_row(colcand, first, all.list[i]);
    if (m_cnt < threadshold) continue;
    genes_order.push_back(i);
    count++;
    candidates[i] = false;
  }
  return count;
}

int add_reverse(const DiscreteArrayListWithSymbols& all, std::vector<int>& genes_reverse, std::vector<char>& candidates,
                const std::list<std::size_t>& colcand, const DiscreteArray& first, double threadshold) {
  std::size_t rows = all.list.size();
  int count = 0;
  for (std::size_t i = 0; i < rows; i++) {
    if (!candidates[i]) continue;
    int m_cnt = internal::reverse_row(colcand, first, all.list[i], all.symbols);
    if (m_cnt < threadshold) continue;
    genes_reverse.push_back(i);
    count++;
    candidates[i] = false;
  }
  return count;
}

/************************************************************************/
std::vector<Block> cluster(const DiscreteArrayListWithSymbols& all, const std::vector<Edge *>& el,
                           std::size_t COL_WIDTH, double TOLERANCE, bool IS_cond, bool IS_area,
                           bool IS_pvalue, std::size_t SCH_BLOCK, int RPT_BLOCK, double FILTER, double f, bool verbose) {
  std::vector<Block> bb;
  std::size_t rows = all.list.size();
  std::size_t cols = all.list[0].size();

  std::set<int> allincluster;
  for (std::vector<Edge *>::const_iterator it = el.begin(); it != el.end(); ++it) {
    const Edge* e = *it;
    /* check if both genes already enumerated in previous blocks */
    bool flag = true;
    /* speed up the program if the rows bigger than 200 */
    if (rows > 250) {
      if (allincluster.find(e->gene_one) != allincluster.end() && allincluster.find(e->gene_two) != allincluster.end()) flag = false;
    } else flag = internal::check_seed(e, bb, rows);
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
    int cand_threshold = static_cast<int>(std::floor(COL_WIDTH * TOLERANCE));
    if (cand_threshold < 2) cand_threshold = 2;
    /* maintain a candidate list to avoid looping through all rows */
    std::vector<char> candidates(rows, true);
    candidates[e->gene_one] = candidates[e->gene_two] = false;
    std::size_t components = 2;
    /* expansion step, generate a bicluster without noise */
    std::vector<long double> pvalues;
    pvalues.reserve(rows);
    internal::block_init(all.list, b, genes_order, scores, candidates, cand_threshold, components, pvalues, IS_cond, COL_WIDTH,
                         IS_area);
    /* track back to find the genes by which we get the best score*/
    std::size_t k;
    for (k = 0; k < components; k++) {
      if (IS_pvalue && (pvalues[k] == b.pvalue && k >= 2 && scores[k] != scores[k + 1])) break;
      if (scores[k] == b.score && (k + 1 == scores.size() || scores[k + 1] != b.score)) break;
    }
    components = k + 1;
    std::fill(candidates.begin(), candidates.end(), true);
    for (std::size_t ki = 0; ki < components - 1; ki++) {
      internal::seed_update(all.list[genes_order[ki]], profile);
      candidates[genes_order[ki]] = false;
    }
    candidates[genes_order[k]] = false;
    genes_order.resize(k + 1);
    std::list<std::size_t> colcand;
    /* add columns satisfy the conservative r */
    int cnt = internal::seed_current_modify(all.list[genes_order[k]], colcand, components, profile, TOLERANCE);

    const DiscreteArray& first = all.list[e->gene_one];
    double threadshold = std::floor(cnt * TOLERANCE);
    /* add some new possible genes */
    components += add_intersect(all, genes_order, candidates, colcand, first, threadshold);
    /* add genes that negative regulated to the consensus */
    components += add_reverse(all, genes_reverse, candidates, colcand, first, threadshold);
    /* store gene arrays inside block */
    internal::scan_block(all.list, all.symbols, genes_order, b, TOLERANCE);
    if (b.block_cols() == 0) continue;
    if (IS_pvalue) b.score = static_cast<int>(-(100 * log(b.pvalue)));
    else b.score = components * b.block_cols();
    for (std::vector<int>::iterator iterator = genes_order.begin(); iterator != genes_order.end(); ++iterator) {
      b.genes_order.insert(*iterator);
      allincluster.insert(*iterator);
    }
    for (std::vector<int>::iterator iterator = genes_reverse.begin(); iterator != genes_reverse.end(); ++iterator) {
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
  return internal::report_blocks(bb, RPT_BLOCK, FILTER);
}

std::vector<Block> read_and_solve_blocks(const DiscreteArrayListWithSymbols& all, double TOLERANCE,
    const std::vector<std::vector<char>>& RowxNumber, const std::vector<std::vector<char>>& NumberxCol) {
  std::size_t rows = all.list.size();

  std::vector<Block> bb;

  for (std::size_t number = 0; number < RowxNumber[0].size(); number++) {
    Block b;
    for (std::size_t row = 0; row < RowxNumber.size(); row++) {
      if (RowxNumber[row][number]) b.genes_order.insert(row);
    }
    for (std::size_t col = 0; col < NumberxCol[number].size(); col++) {
      if (NumberxCol[number][col]) b.conds.insert(col);
    }

    {
      std::vector<int> genes_order, genes_reverse;
      genes_order.reserve(rows);
      genes_reverse.reserve(rows);
      /* maintain a candidate list to avoid looping through all rows */
      std::vector<char> candidates(rows, true);
      for (std::size_t row = 0; row < RowxNumber.size(); row++) {
        if (RowxNumber[row][number]) {
          candidates[row] = false;
          genes_order.push_back(row);
        }
      }

      std::list<std::size_t> colcand;
      for (std::size_t col = 0; col < NumberxCol[number].size(); col++) {
        if (NumberxCol[number][col]) {
          colcand.push_back(col);
        }
      }

      double threadshold = std::floor(colcand.size() * TOLERANCE);

      const DiscreteArray& first = all.list[*b.genes_order.begin()];

      /* add some new possible genes */
      add_intersect(all, genes_order, candidates, colcand, first, threadshold);
      /* add genes that negative regulated to the consensus */
      add_reverse(all, genes_reverse, candidates, colcand, first, threadshold);

      for (auto it = genes_order.begin(); it != genes_order.end(); ++it) {
        b.genes_order.insert(*it);
      }
      for (auto it = genes_reverse.begin(); it != genes_reverse.end(); ++it) {
        b.genes_reverse.insert(*it);
      }

      for (auto it = colcand.begin(); it != colcand.end(); ++it) {
        b.conds.insert(*it);
      }
    }

    bb.push_back(b);
  }

  return bb;// internal::report_blocks(bb, RPT_BLOCK, FILTER);
}
