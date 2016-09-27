#ifndef __TENSORQUADRATURE__CLASS__
#define __TENSORQUADRATURE__CLASS__

#include <Eigen/Dense>

#include "UnivariateQuadrature.hpp"

class TensorProductQuadrature {
 public:
  TensorProductQuadrature(void);
  TensorProductQuadrature(const Eigen::VectorXi &lvl, UnivariateQuadrature &Q);
  void initQuadrature(const Eigen::VectorXi &lvl, UnivariateQuadrature &Q);
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
