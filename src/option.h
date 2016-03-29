#ifndef OPTION_H
#define OPTION_H

struct Option {
  Option(bool _p, bool _s, bool _c, bool _f) : p(_p), s(_s), c(_c), filter_1xn_nx1(_f){}
  Option() : p(false), s(false), c(false), filter_1xn_nx1(false) {}

  bool p;
  bool s;
  bool c;
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
