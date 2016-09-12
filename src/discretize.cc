#include "discretize.h"

#include <cmath>
#include <algorithm>

continuous quantile_from_sorted_data(const std::vector<continuous> &sorted_data, const std::size_t n, const double f) {
  /*std::floor function returns the largest integral value less than or equal to x*/
  int i = static_cast <int>(std::floor((n - 1) * f));
  continuous delta = static_cast<continuous>((n - 1) * f - i);
  return (1 - delta) * sorted_data[i] + delta * sorted_data[i + 1];
}

short dis_value(const float current, const short divided, const std::vector<continuous> &small, const int cntl,
  const std::vector<continuous> &big, const int cntu) {
  continuous d_space = static_cast<continuous>(1.0 / divided);
  for (short i = 0; i < divided; i++) {
    if (cntl > 0 && current <= quantile_from_sorted_data(small, cntl, static_cast<continuous>(d_space * (i + 1))))
      return -i - 1;
    if (cntu > 0 && current >= quantile_from_sorted_data(big, cntu, static_cast<continuous>(1.0 - d_space * (i + 1))))
      return i + 1;
  }
  return 0;
}

rule do_work(const std::vector<std::vector<continuous>>& arr, const double f, const short divided,
  DiscreteArrayList & arr_d, const std::size_t row, std::vector<continuous> & rowdata, std::vector<continuous> & big, std::vector<continuous> & small)
{
  for (std::size_t col = 0; col < arr[0].size(); col++)
    rowdata[col] = arr[row][col];
  std::sort(rowdata.begin(), rowdata.end());
  float f1 = quantile_from_sorted_data(rowdata, arr[0].size(), 1 - f);
  float f2 = quantile_from_sorted_data(rowdata, arr[0].size(), f);
  float f3 = quantile_from_sorted_data(rowdata, arr[0].size(), 0.5);

  float upper, lower;

  if (f1 - f3 >= f3 - f2) {
    upper = 2 * f3 - f2;
    lower = f2;
  } else {
    upper = f1;
    lower = 2 * f3 - f1;
  }
  std::size_t cntu = 0;
  std::size_t cntl = 0;
  for (std::size_t i = 0; i < arr[0].size(); i++) {
    if (rowdata[i] < lower) {
      small[cntl] = rowdata[i];
      cntl++;
    }
    if (rowdata[i] > upper) {
      big[cntu] = rowdata[i];
      cntu++;
    }
  }
  for (std::size_t col = 0; col < arr[0].size(); col++)
    arr_d[row][col] = dis_value(arr[row][col], divided, small, cntl, big, cntu);

  rule rule;
  rule.lower = lower;
  rule.upper = upper;
  rule.cntl = cntl;
  rule.cntu = cntu;
  return rule;
}

DiscreteArrayList discretize(const std::vector<std::vector<continuous>> &arr,
  const short divided, const double f) {
  DiscreteArrayList arr_d(arr.size(), DiscreteArray(arr[0].size()));
  std::vector<continuous> rowdata(arr[0].size());
  std::vector<continuous> big(arr[0].size()), small(arr[0].size());
  for (std::size_t row = 0; row < arr.size(); row++) {
    do_work(arr, f, divided, arr_d, row, rowdata, big, small);
  }
  return arr_d;
}

DiscreteArrayList discretize(const std::vector<std::vector<continuous>> &arr,
  const short divided, const double f, std::vector<rule> &genes_rules) {
  DiscreteArrayList arr_d(arr.size(), DiscreteArray(arr[0].size()));
  std::vector<continuous> rowdata(arr[0].size());
  std::vector<continuous> big(arr[0].size()), small(arr[0].size());
  for (std::size_t row = 0; row < arr.size(); row++) {
    genes_rules.push_back(do_work(arr, f, divided, arr_d, row, rowdata, big, small));
  }
  return arr_d;
}
