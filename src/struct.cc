#include "struct.h"

/* Return the number of common components between two arrays */
int count_intersect(const std::set<int> & ds1, const std::set<int> & ds2) {
  int cnt = 0;
  std::set<int>::const_iterator first1 = ds1.begin();
  std::set<int>::const_iterator last1 = ds1.end();
  std::set<int>::const_iterator first2 = ds2.begin();
  std::set<int>::const_iterator last2 = ds2.end();
  while (first1 != last1 && first2 != last2) {
    if (*first1 < *first2) ++first1;
    else if (*first2 < *first1) ++first2;
    else {
      ++cnt;
      ++first1;
      ++first2;
    }
  }
  return cnt;
}
