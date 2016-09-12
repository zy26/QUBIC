#include "charset.h"

#include <climits> // SHRT_MAX
#include <cstdio>
#include <cstddef> // std::size_t

short charset_add(DiscreteArray &ar, const short &s, short *bb) {
  /*A signed short can hold all the values between SHRT_MIN  and SHRT_MAX inclusive.SHRT_MIN is required to be -32767 or less,SHRT_MAX must be at least 32767*/
  int ps = s + SHRT_MAX;
  if (bb[ps] < 0) {
    bb[ps] = static_cast<short>(ar.size());
    ar.push_back(s);
  }
  return bb[ps];
}

DiscreteArrayListWithSymbols make_charsets_d(const DiscreteArrayList &arr, bool verbose) {
  DiscreteArrayListWithSymbols all;
  all.list.resize(arr.size(), DiscreteArray(arr[0].size()));
  short bb[USHRT_MAX];
  std::fill_n(bb, USHRT_MAX, -1);
  charset_add(all.symbols, 0, bb);
  for (std::size_t i = 0; i < arr.size(); i++)
    for (std::size_t j = 0; j < arr[0].size(); j++)
      all.list[i][j] = charset_add(all.symbols, arr[i][j], bb);
  if (verbose) fprintf(stdout, "Discretized data contains %d classes with charset [ ", static_cast<unsigned int>(all.symbols.size()));
  for (std::size_t i = 0; i < all.symbols.size(); i++)
    if (verbose) fprintf(stdout, "%d ", all.symbols[i]);
  if (verbose) fprintf(stdout, "]\n");
  return all;
}
