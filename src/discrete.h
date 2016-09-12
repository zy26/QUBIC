#ifndef DISCRETE_H
#define DISCRETE_H

#include <vector>

typedef std::vector<short> DiscreteArray;
typedef std::vector<DiscreteArray> DiscreteArrayList;

typedef std::vector<short> Symbols;

class DiscreteArrayListWithSymbols {
public:
  DiscreteArrayList list;
  Symbols symbols;
private:
};
#endif // DISCRETE_H
