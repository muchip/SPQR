#ifndef __HCINDEXSET__CLASS__
#define __HCINDEXSET__CLASS__

#include "SparseIndexSet.hpp"

#include <Eigen/Dense>

/** \brief Hyperbolic cross is implemented as a simple addon for the
*          sparseIndexSet class which provides the corresponding comparison
*          functor. Weights my be unsorted in HC case and no sorting is
*          performed.
*/
class HCindexSet : public SparseIndexSet {
 public:
  void computeIndexSet(int q, const Eigen::VectorXd &w) {
    SparseIndexSet::computeIndexSet(q, w.size(), HCindexSet::cpFun(w));
  };

 protected:
  struct cpFun {
    Eigen::VectorXd _w;
    cpFun(const Eigen::VectorXd &w) : _w(w){};
    double operator()(const Eigen::VectorXi &alpha) const {
      double skap = 1;
      for (int i = 0; i < _w.size(); ++i)
        skap *= std::pow(alpha(i) + 1., _w(i));
      return skap - 1.;
    }
  };
};
#endif
