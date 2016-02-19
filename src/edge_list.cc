#include "edge_list.h"
#include "config.h"
#include "fib.h"
#include <cassert> // assert
#include <cstdio>
#include <algorithm>


int str_intersect_r(const DiscreteArray &s1, const DiscreteArray &s2) {
  assert(s1.size() == s2.size());
  int common_cnt = 0;
  /* s1 and s2 of equal length, so we check s1 only */
  for (size_t i = 0; i < s1.size(); i++)
    if ((s1[i] != 0) && (s1[i] == s2[i])) // Changed order by zy26
      common_cnt++;
  return common_cnt;
}

static int edge_cmpr(void *a, void *b) {
  int score_a, score_b;
  score_a = ((Edge *)a)->score;
  score_b = ((Edge *)b)->score;
  if (score_a < score_b) return -1;
  if (score_a == score_b) return 0;
  return 1;
}

static int fh_insert_fixed(fibheap *h, Edge *data, int cur_min) {
  if (h->fh_n < HEAP_SIZE)
    fh_insert(h, data);
  else {
    if (data->score > cur_min) {
      /* Remove least value and renew */
      fh_extractmin(h);
      fh_insert(h, data);
      /* Keep a memory of the current min */
      cur_min = ((Edge *)fh_min(h))->score;
    }
  }
  return cur_min;
}

static void fh_dump(fibheap *h, std::vector<Edge *> &data_array, int min_score) {
  int i;
  int n = h->fh_n;
  for (i = n - 1; i >= 0; i--) {
    Edge *t = (Edge *)fh_min(h);
    if (t->score > min_score) break;
    fh_extractmin(h);
  }
  data_array.resize(i + 1);
  for (; i >= 0; i--)
    data_array[i] = (Edge *)fh_extractmin(h);
}

int EdgeList::get_key(const Edge *s) {
  return s->score - col_width;
}

int get_key1(const Edge &s) {
  return s.score;
}

struct CompEventByPtr {
  bool operator()(const Edge *pEvent1, const Edge *pEvent2) const {
    return pEvent1->score >= pEvent2->score;
  }
};

const std::vector<Edge *> &EdgeList::get_edge_list() const {
  return edge_list;
}
EdgeList::EdgeList(const DiscreteArrayList &arr_c, size_t &COL_WIDTH, bool verbose) {
  Edge *edge;
  int cnt;
  /* Allocating heap structure */
  fibheap *heap = fh_makeheap();
  fh_setcmp(heap, edge_cmpr);
  /* Generating seed list and push into heap */
  if (verbose) fprintf(stdout, "Generating seed list (minimum weight %d)\n", static_cast<unsigned int>(COL_WIDTH));
  int min_score = COL_WIDTH - 1;
  /* iterate over all genes to retrieve all edges */
  for (size_t i = 0; i < arr_c.size(); i++)
    for (size_t j = i + 1; j < arr_c.size(); j++) {
      cnt = str_intersect_r(arr_c[i], arr_c[j]);
      if (cnt <= min_score) continue;
      edge = new Edge();
      edge->gene_one = i;
      edge->gene_two = j;
      edge->score = cnt;
      min_score = fh_insert_fixed(heap, edge, min_score);
    }
  if (heap->fh_n == 0) {
    fprintf(stderr, "[Error] Not enough overlap between genes");
    throw 1.0;
  }
  /* sort the seeds */
  if (verbose) fprintf(stdout, "%d seeds generated\n", heap->fh_n);
  fh_dump(heap, edge_list, min_score);
  if (verbose) fprintf(stdout, "%d seeds dumped\n", static_cast<unsigned int>(edge_list.size()));
}

EdgeList::~EdgeList() {
  for (size_t i = 0; i < edge_list.size(); i++)
    delete(edge_list[i]);
}
