#ifndef DISCRETE_H
#define DISCRETE_H

#include <vector>

typedef short discrete;
typedef std::vector<discrete> DiscreteArray;
typedef std::vector<DiscreteArray> DiscreteArrayList;

typedef std::vector<discrete> Symbols;

class DiscreteArrayListWithSymbols
{
public:
  DiscreteArrayList list;
  Symbols symbols;
private:

};
#endif // DISCRETE_H
