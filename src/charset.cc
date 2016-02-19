#include "charset.h"

#include <climits> // SHRT_MAX
#include <cstdio>
#include <cstddef> // size_t

discrete charset_add(std::vector<discrete> &ar, const discrete &s, discrete *bb) {
  /*A signed short can hold all the values between SHRT_MIN  and SHRT_MAX inclusive.SHRT_MIN is required to be -32767 or less,SHRT_MAX must be at least 32767*/
  int ps = s + SHRT_MAX;
  if (bb[ps] < 0) {
    bb[ps] = static_cast<discrete>(ar.size());
    ar.push_back(s);
  }
  return bb[ps];
}

DiscreteArrayListWithSymbols make_charsets_d(const std::vector<std::vector<discrete>> &arr, bool verbose) {
  DiscreteArrayListWithSymbols all;
  all.list.resize(arr.size(), DiscreteArray(arr[0].size()));
  discrete bb[USHRT_MAX];
  std::fill_n(bb, USHRT_MAX, -1);
  charset_add(all.symbols, 0, bb);
  for (size_t i = 0; i < arr.size(); i++)
    for (size_t j = 0; j < arr[0].size(); j++)
      all.list[i][j] = charset_add(all.symbols, arr[i][j], bb);
  if (verbose) fprintf(stdout, "Discretized data contains %d classes with charset [ ", static_cast<unsigned int>(all.symbols.size()));
  for (size_t i = 0; i < all.symbols.size(); i++)
    if (verbose) fprintf(stdout, "%d ", all.symbols[i]);
  if (verbose) fprintf(stdout, "]\n");
  return all;
}
