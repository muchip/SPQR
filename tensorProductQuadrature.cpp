#include "tensorProductQuadrature.hpp"

/**    \brief default constructor without arguments
*
*/
tensorProductQuadrature::tensorProductQuadrature(void) {}

/**    \brief constructor with quadrature degrees calls init_quadrature
*      \param[in] degs maximum quadrature degrees in each dimension
*
*/
tensorProductQuadrature::tensorProductQuadrature(const Eigen::VectorXi &lvl,
                                                 univariateQuadrature &Q) {
  init_quadrature(lvl, Q);
}

/**    \brief initializes the tensor product quadrature based on the
*             Gauss Legendre quadrature (C routine)
*      \param[in] degs maximum quadrature degrees in each dimension
*
*/
void tensorProductQuadrature::init_quadrature(const Eigen::VectorXi &lvl,
                                              univariateQuadrature &Q) {
  Eigen::VectorXi base;
  Eigen::VectorXi nPtsInQuad;
  int remainder = 0;
  int quotient = 0;
  int maxLvl = 0;

  for (int i = 0; i < lvl.size(); ++i)
    if (maxLvl < lvl(i)) maxLvl = lvl(i);

  Q.resizeQuadrature(maxLvl);
  nPtsInQuad.resize(lvl.size());
  const std::vector<Quadrature> &refQ = Q.get_Q();

  for (int i = 0; i < lvl.size(); ++i) nPtsInQuad(i) = refQ[lvl(i)].w.size();

  _nPts = nPtsInQuad.prod();
  _points.resize(lvl.size(), _nPts);
  _weights.resize(_nPts);

  base = computeBase(nPtsInQuad);
  for (int i = 0; i < _nPts; ++i) {
    _weights(i) = 1.;
    remainder = i;
    for (int j = 0; j < lvl.size(); ++j) {
      quotient = remainder / base(j);
      _points(j, i) = refQ[lvl(j)].xi(quotient);
      _weights(i) *= refQ[lvl(j)].w(quotient);
      remainder -= base(j) * quotient;
    }
  }
}

/**   \brief this function needs a comment
*
*/
Eigen::VectorXi tensorProductQuadrature::computeBase(
    const Eigen::VectorXi &lvl) {
  Eigen::VectorXi base;
  base = Eigen::VectorXi::Ones(lvl.size());
  for (int i = (int)lvl.size() - 1; i >= 0; --i)
    base(lvl.size() - i - 1) = lvl.tail(i).prod();
  return base;
}

/**   \brief make an educated guess, what this function does...
*
*/
const Eigen::VectorXd &tensorProductQuadrature::get_weights(void) const {
  return _weights;
}

/**   \brief make an educated guess, what this function does...
*
*/
const Eigen::MatrixXd &tensorProductQuadrature::get_points(void) const {
  return _points;
}

/**   \brief make an educated guess, what this function does...
*
*/
int tensorProductQuadrature::get_nPts(void) const { return _nPts; }

