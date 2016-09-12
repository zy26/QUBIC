#ifndef EDGE_LIST_H
#define EDGE_LIST_H

#include "discrete.h"
#include <limits>
#include <cstddef> // size_t
#include <cassert>
#include <algorithm>

/* edge between two genes */
struct Edge {
  unsigned int gene_one;
  unsigned int gene_two;
  int score;
};

template <typename T>
class AdjMatrix {
  std::vector<T> matrix_;
  unsigned int size_;

public:
  AdjMatrix(unsigned int size) : size_(size) {
    matrix_.resize(size_* (size_ + 1) / 2);
  }

  std::size_t get_index(unsigned x, unsigned y) const {
    return (y * (2 * size_ - y + 1)) / 2 + (x - y);
  }

  T& operator()(unsigned int x, unsigned int y) {
    assert(x < size_ && y < size_);
    if (x >= y) return matrix_[get_index(x, y)];
    return matrix_[get_index(y, x)];
  }
  const T& operator()(unsigned int x, unsigned int y) const {
    assert(x < size_ && y < size_);
    if (x >= y) return matrix_[get_index(x, y)];
    return matrix_[get_index(y, x)];
  }
};

inline unsigned str_intersect_r(const DiscreteArray &s1, const DiscreteArray &s2) {
  assert(s1.size() == s2.size());
  int common_cnt = 0;
  /* s1 and s2 of equal length, so we check s1 only */
  for (std::size_t i = 0; i < s1.size(); i++)
    if (s1[i] != 0 && s1[i] == s2[i]) // Changed order by zy26
      common_cnt++;
  return common_cnt;
}

class CountHelper {
protected:
  const DiscreteArrayList &arr_;
private:
  const std::size_t col_width_;
public:
  virtual std::size_t col_width() const {
    return col_width_;
  }
  virtual ~CountHelper() {
  }

  explicit CountHelper(const DiscreteArrayList& arr_c, std::size_t col_width) : arr_(arr_c), col_width_(col_width) { }

  std::size_t size() const {
    return arr_.size();
  }
  // I want to change following function form virtual to abstract (=0), but unfortunately it will slow down the program under ms windows (other platform untested).
  // I don't known why.
  virtual int get_weight(std::size_t i, std::size_t j) const {
    return str_intersect_r(arr_[i], arr_[j]);
  }

  virtual void Update(std::vector<Edge*>& /*edges*/) const {};
};

class CountHelperRealTime : public CountHelper {
  explicit CountHelperRealTime(const DiscreteArrayList& arr_c, std::size_t col_width) : CountHelper(arr_c, col_width) { }
};

class CountHelperSaved : public CountHelper {
protected:
  std::vector<unsigned> intersects_;
public:
  virtual ~CountHelperSaved() {
  }

  explicit CountHelperSaved(const DiscreteArrayList& arr_c, std::size_t col_width) : CountHelper(arr_c, col_width), intersects_(arr_c.size() * (arr_c.size() - 1) / 2) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; static_cast<std::size_t>(i) < arr_c.size(); i++)
      for (std::size_t j = i + 1; j < arr_c.size(); j++)
        intersects_[j * (j - 1) / 2 + i] = str_intersect_r(arr_c[i], arr_c[j]); // if not compress, it will save some time, but memory cannot be alloc sometimes.
  }

  int get_weight(std::size_t i, std::size_t j) const override {
    assert(i < j);
    return intersects_[j * (j - 1) / 2 + i];
  }
};

class CountHelperRanked : public CountHelperSaved {
  struct mycomparison {
    bool operator() (unsigned* lhs, unsigned* rhs) const {
      return *lhs < *rhs;
    }
  };
  std::size_t col_width_;
protected:
  std::vector<unsigned> weight_;
public:
  std::size_t col_width() const override {
    return col_width_;
  }

  virtual ~CountHelperRanked() {
  }

  explicit CountHelperRanked(const DiscreteArrayList& arr_c, std::size_t col_width) : CountHelperSaved(arr_c, col_width), weight_(intersects_) {
    assert(weight_.size() <= std::numeric_limits<unsigned int>::max());
    unsigned int size = static_cast<unsigned int>(weight_.size());
    std::vector<unsigned*> pintArray((size));
    for (unsigned int i = 0; i < size; ++i) {
      pintArray[i] = &weight_[i];
    }

    std::sort(pintArray.begin(), pintArray.end(), mycomparison());

    // Dereference the pointers and assign their sorted position. not deal tie
    for (unsigned int i = 0; i < size; ++i) {
      *pintArray[i] = i + 1;
    }
    col_width_ = col_width;
  }

  void Update(std::vector<Edge*>& edges) const override {
    for (auto it = edges.begin(); it != edges.end(); ++it) {
      int i = (*it)->gene_one;
      int j = (*it)->gene_two;
      (*it)->score = weight_[j * (j - 1) / 2 + i];
    }
  };
};

class WeightedCountHelper : public CountHelperRanked {
  const std::vector<std::vector<float>>& input_weights_;
public:
  explicit WeightedCountHelper(const DiscreteArrayList& arr_c, const std::vector<std::vector<float>>& weights, std::size_t col_width) : CountHelperRanked(arr_c, col_width), input_weights_(weights) {}

  void Update(std::vector<Edge*>& edges) const override {
    for (auto it = edges.begin(); it != edges.end(); ++it) {
      int i = (*it)->gene_one;
      int j = (*it)->gene_two;
      (*it)->score = weight_[j * (j - 1) / 2 + i] + static_cast<int>(input_weights_[i][j]);;
    }
  };
};

class EdgeList {
  std::vector<Edge *> edge_list_;
public:
  const std::vector<Edge *> &get_seeds() const;
  EdgeList(const CountHelper& countHelper, bool verbose);
  ~EdgeList();
};

#endif
