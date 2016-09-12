#include "edge_list.h"
#include "config.h"
#include "fib.h"
#include <cassert> // assert
#include <cstdio>
#include <algorithm>

static int edge_cmpr(void *a, void *b) {
  int score_a = static_cast<Edge *>(a)->score;
  int score_b = static_cast<Edge *>(b)->score;
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
      cur_min = static_cast<Edge *>(fh_min(h))->score;
    }
  }
  return cur_min;
}

static void fh_dump(fibheap *h, std::vector<Edge *> &data_array, int min_score) {
  int i;
  int n = h->fh_n;
  for (i = n - 1; i >= 0; i--) {
    Edge *t = static_cast<Edge *>(fh_min(h));
    if (t->score > min_score) break;
    fh_extractmin(h);
  }
  data_array.resize(i + 1);
  for (; i >= 0; i--)
    data_array[i] = static_cast<Edge *>(fh_extractmin(h));
}

struct CompEventByPtr {
  bool operator()(const Edge *pEvent1, const Edge *pEvent2) const {
    return pEvent1->score > pEvent2->score;
  }
};

const std::vector<Edge *> &EdgeList::get_seeds() const {
  return edge_list_;
}

EdgeList::EdgeList(const CountHelper& countHelper, bool verbose) {
  std::size_t col_width = countHelper.col_width();
  Edge *edge;
  int cnt;
  /* Allocating heap structure */
  fibheap *heap = fh_makeheap();
  fh_setcmp(heap, edge_cmpr);
  /* Generating seed list and push into heap */
  if (verbose) fprintf(stdout, "Generating seed list (minimum weight %d)\n", static_cast<unsigned int>(col_width));
  int min_score = col_width - 1;
  /* iterate over all genes to retrieve all edges */
  for (std::size_t i = 0; i < countHelper.size(); i++)
    for (std::size_t j = i + 1; j < countHelper.size(); j++) {
      cnt = countHelper.get_weight(i, j);
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
  fh_dump(heap, edge_list_, min_score);
  countHelper.Update(edge_list_);
  std::stable_sort(edge_list_.begin(), edge_list_.end(), CompEventByPtr()); // I guess we can remove stable sort or fib someday
  if (verbose) fprintf(stdout, "%d seeds dumped\n", static_cast<unsigned int>(edge_list_.size()));
}

EdgeList::~EdgeList() {
  for (std::size_t i = 0; i < edge_list_.size(); i++)
    delete(edge_list_[i]);
}
