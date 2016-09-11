#include "TDindexSet.hpp"

/**    \brief copies the vector gamma to the class pendant and sorts it via
*             sortGamma
*      \param[in] gamma anisotropies in each dimension
*
*/
void TDindexSet::set_w(const Eigen::VectorXd &w) {
  _w = w;
  TDindexSet::init_sortW();
  // reorder _w according to sortW;
  for (int i = 0; i < w.size(); ++i) _w(i) = w(_sortW(i));
}

/**    \brief initializes the sparse anisotropic index set
*      \param[in] q is the maximum allowed value of the weighted sum of
*                 gamma times index
*
*/
void TDindexSet::comp_indexSet(int q, const Eigen::VectorXd &w) {
  int k = 0;
  Eigen::VectorXi currInd(w.size());

  // set max level for the sparse index set
  _q = q;
  // set dimension of the sparse grid
  _dim = (int)w.size();
  // set weight vector
  TDindexSet::set_w(w);

  _alpha.resize(_dim, __MEMCHUNKSIZE__);
  _cw.resize(__MEMCHUNKSIZE__);
  _alpha.setZero();
  _cw.setZero();
  _myOnes = Eigen::VectorXi::Ones(_dim);

  currInd.setZero();
  k = 0;
  // test if 0 is in Xw, else index set empty due to downward closedness
  if (0 <= _q) {
    _sumW = _w.sum();
    if (_sumW > _q) _cw(0) = TDindexSet::cw_alpha(_q, 0, 1, 1);
    if (_cw(0)) ++k;
    TDindexSet::Yw_alpha(0, &k, _q, currInd);
  }

  _alpha.conservativeResize(_dim, k);
  _cw.conservativeResize(k);
}

/**    \brief sorts the values in _w with inreasing magnitude
*      \param[out] sortW permutation vector
*
*/
void TDindexSet::init_sortW(void) {
  _sortW = Eigen::ArrayXi::LinSpaced((int)_w.size(), 0, (int)_w.size() - 1);

  std::sort(&(_sortW(0)), &(_sortW(0)) + _sortW.size(),
            TDindexSet::myCompareInc<Eigen::VectorXd>(_w));
}

/**    \brief computes indices in weighted sparse grid space (recursive)
*
*/
void TDindexSet::Yw_alpha(int maxBit, int *k, double q,
                          Eigen::VectorXi &currInd) {
  int cw = 0;
  double scap = 0;
  for (int i = maxBit; i < _dim; ++i) {
    ++currInd(i);
    q -= _w(i);
    if (q >= 0) {
      scap = _w.dot(currInd.cast<double>());
      if (scap > _q - _sumW)
        cw = TDindexSet::cw_alpha(_q - scap, 0, 1, 1);
      else
        cw = 0;
      if (cw) {
        if (_alpha.cols() <= *k) {
          _alpha.conservativeResize(_dim, _alpha.cols() + __MEMCHUNKSIZE__);
          _cw.conservativeResize(_cw.size() + __MEMCHUNKSIZE__);
        }
        _alpha.col(*k) = currInd;
        _cw(*k) = cw;
        ++(*k);
      }
      Yw_alpha(i, k, q, currInd);
      q += _w(i);
      --currInd(i);
    } else {  // this is the major difference to the base class, we may break
              // here due to the increasingly ordered weights
      q += _w(i);
      --currInd(i);
      break;
    }
  }
}

/**    \brief computes weights for the tensor product quadrature (recursive)
*
*/
int TDindexSet::cw_alpha(double q, int maxBit, int cw, int lvl) {
  for (int i = maxBit; i < _dim; ++i) {
    q -= _w(i);
    if (q >= 0) {
      if (lvl % 2)
        --cw;
      else
        ++cw;
      cw = TDindexSet::cw_alpha(q, i + 1, cw, lvl + 1);
      q += _w(i);
    } else {  // this is the major difference to the base class, we may break
              // here due to the increasingly ordered weights
      q += _w(i);
      break;
    }
  }

  return cw;
}

/**   \brief make an educated guess, what this function does...
*
*/
const Eigen::VectorXi &TDindexSet::get_sortW(void) const { return _sortW; }
