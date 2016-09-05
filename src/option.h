#ifndef OPTION_H
#define OPTION_H

struct Option {
  Option(bool _p, bool _s, bool _c, bool _f) : pvalue_(_p), area_(_s), cond_(_c), filter_1xn_nx1(_f) {}
  Option() : pvalue_(false), area_(false), cond_(false), filter_1xn_nx1(false) {}

  bool pvalue_;
  bool area_;
  bool cond_;
  bool filter_1xn_nx1;
};

struct DOption {
  double c;
  int o;
  double f;
  int k;
  Option option;
  bool verbose;
};

#endif
