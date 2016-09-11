#ifndef __TENSORQUADRATURE__CLASS__
#define __TENSORQUADRATURE__CLASS__

#include <Eigen/Dense>

#include "univariateQuadrature.hpp"

class tensorProductQuadrature {
 public:
  tensorProductQuadrature(void);
  tensorProductQuadrature(const Eigen::VectorXi &lvl, univariateQuadrature &Q);
  void init_quadrature(const Eigen::VectorXi &lvl, univariateQuadrature &Q);
  const Eigen::VectorXd &get_weights(void) const;
  const Eigen::MatrixXd &get_points(void) const;
  int get_nPts(void) const;
  Eigen::VectorXi computeBase(const Eigen::VectorXi &lvl);

 protected:
  Eigen::VectorXd _weights;
  Eigen::MatrixXd _points;
  int _nPts;
};

#endif
