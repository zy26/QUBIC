#ifndef DISCRETIZE_H
#define DISCRETIZE_H

#include <vector>
#include <cstddef> // size_t

#include "discrete.h"
#include "struct.h"

struct rule {
  float lower;
  float upper;
  std::size_t cntl;
  std::size_t cntu;
};

DiscreteArrayList discretize(const std::vector<std::vector<continuous>> &arr, const short divided, const double f,
  std::vector<rule> &genes_rules);
DiscreteArrayList discretize(const std::vector<std::vector<continuous>> &arr, const short divided, const double f);
#endif
