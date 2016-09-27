#ifndef __TDINDEXSET__CLASS__
#define __TDINDEXSET__CLASS__

// include for C++ sort function
#include <algorithm>

// Eigen Library
#include <Eigen/Dense>

// include base class
#include "SparseIndexSet.hpp"

#include "CONSTANTS.hpp"


class TDindexSet : public SparseIndexSet {
 public:
  // override base class methods
  void computeIndexSet(int q, const Eigen::VectorXd &w);

  const Eigen::VectorXi &get_sortW(void) const;

 protected:
  // override base class methods
  void combiIndexSet(int maxBit, int *k, double q, Eigen::VectorXi &currInd);
  int combiWeights(double q, int maxBit, int cw, int lvl);

  // new methods and variables related to weight vector
  void init_sortW(void);
  void set_w(const Eigen::VectorXd &w);

  Eigen::VectorXd _w;
  Eigen::VectorXi _sortW;
  double _sumW;
  // comparison functors for C++ sort routine
  template <typename T>
  struct myCompareInc {
    const T &m;
    myCompareInc(T &p) : m(p){};

    bool operator()(const int &i, const int &j) { return (m(i) < m(j)); }
  };
};
#endif

#if 0
template <typename T>
struct myCompareDec {
    const T &m;
    myCompareDec(T &p) : m(p){};
    
    bool operator()(const int &i, const int &j) { return (m(i) > m(j)); }
};
#endif
