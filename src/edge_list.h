#ifndef EDGE_LIST_H
#define EDGE_LIST_H

#include <cstddef> // size_t

#include "discrete.h"
#include "struct.h"

/* edge between two genes */
typedef struct Edge {
  size_t gene_one;
  size_t gene_two;
  int score;
} Edge;

class EdgeList {
private:
  std::vector<Edge *> edge_list; int col_width;
  int get_key(const Edge* s);
public:
  const std::vector<Edge *> &get_edge_list() const;
  EdgeList(const DiscreteArrayList &, size_t&, bool verbose);
  ~EdgeList();
};

#endif
