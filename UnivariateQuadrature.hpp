#ifndef __UNIVARIATEQUADRATURE__CLASS__
#define __UNIVARIATEQUADRATURE__CLASS__

#include <vector>

#include <Eigen/Dense>

// struct that contains a particular quadrature rule
struct Quadrature {
  Eigen::VectorXd xi;
  Eigen::VectorXd w;
};

/**   \brief basis class for univariate quadrature rules
*
*/
class UnivariateQuadrature {
 public:
  UnivariateQuadrature(void) : _maxLvl(-1){};
  // setter...
  virtual void initQuadrature(int maxLvl) = 0;
  virtual void resizeQuadrature(int maxLvl) = 0;
  // tester...
  virtual void testQuadrature(int maxLvl) = 0;
  // getter...
  int get_maxLvl(void) const { return _maxLvl; };
  const std::vector<Quadrature> &get_Q(void) const {
    return (const std::vector<Quadrature> &)_Q;
  };

 protected:
  int _maxLvl;
  std::vector<Quadrature> _Q;
};

#endif
