#ifndef __SPARSEINDEXSET__CLASS__
#define __SPARSEINDEXSET__CLASS__

// Eigen Library
#include <Eigen/Dense>

#include "CONSTANTS.hpp"

class SparseIndexSet {
 public:
  /**     \brief initializes the sparse anisotropic index set
  *       \param[in] q is the maximum level for the sparse index set
  *       \param[in] dim is the dimension of the sparse index set
  *       \param[in] cpFun is the comparison Function. sparse index set is
  *                  constructed such that cpFun(alpha) <= q for all alpha in
  *                  the sparse index set
  *
  * 	    Additionally, the set of multiindices _alpha and the weights _cw
  *     	 are provided globally.
  *
  * 	    The function comp_indexSet computes a sparse multi index set
  * 	    _alpha\subset\N^n by the function call
  *       sparseIndexSet::Yw_alpha(0, &n, cpFun)
  * 	    and corresponding weights _cw \in \Z^n to each multi index.
  * 	    To determine the weights, we have to check for each multi index
  *       alpha contained in _alpha if alpha + beta is in _alpha for each beta 
  *       in {0,1}^{n}. If this is the case, we update the weight by 
  *       (-1)^|beta|. Of course, the condition is always fulfilled for the 
  *       multiindex beta = 0. Therefore, we skip this multiindex and initialize
  *       the weight cw with 1 which is associated with the second variable in
  *       the function call to cw_alpha.
  * 	    For a multi index alpha, the corresponding weight has only to be
  *       computed if the multiindex alpha+1 is not contained in _alpha.
  *       Otherwise, the corresponding weight is 0. With the index set, the
  *       weights and the unidirectional quadrature formulas at hand, the
  *       integral I(f) is approximated by the combination technique formula
  * 	    I(f) \approx \sum_{i =0}^{_alpha.cols-1} _cw(i)*Q_{alpha.col(i)}.
  * 	    Herein, _alpha.cols determines the number of multiindices in _alpha
  *       and alpha.col(i) gives the i-th index in _alpha. The tensor product
  *       quadrature formula Q_{alpha.col(i)} is specified by the unidirectional
  *       quadrature formulas.
  */
  template <class Compare>
  void computeIndexSet(int q, int dim, const Compare &cpFun) {
    int k = 0;
    Eigen::VectorXi currInd(dim);
    // set max level for the sparse index set
    _q = q;
    // set dimension of the sparse grid
    _dim = dim;
    // allocate memory for _alpha and _cw to keep the overhead at a low level
    // blocks of size __MEMCHUNKSIZE__ are allocated at once
    _alpha.resize(_dim, __MEMCHUNKSIZE__);
    _cw.resize(__MEMCHUNKSIZE__);
    // set memory to 0
    _alpha.setZero();
    _cw.setZero();
    // get a vector with all ones to reduce the overhead for the constructor of
    // Eigen::VectorXi::Ones(_dim);
    _myOnes = Eigen::VectorXi::Ones(_dim);
    // start from multi index 0\in Xw_alpha
    currInd.setZero();
    k = 0;
    // test if 0 is in Xw_alpha, else index set empty due to downward closedness
    if (cpFun(currInd) <= _q) {
      // check if _cw(0)\neq 0
      if (cpFun(_myOnes) > _q) _cw(0) = combiWeights(0, 1, 1, currInd, cpFun);
      if (_cw(0)) ++k;
      // compute all other indices in Yw_alpha recursively
      combiIndexSet(0, &k, cpFun, currInd);
    }
    // crop memory for _alpha and _cw to actual size
    _alpha.conservativeResize(_dim, k);
    _cw.conservativeResize(k);
  }

  /**   \brief make an educated guess, what this function does...
  *
  */
  const Eigen::MatrixXi &get_alpha(void) const { return _alpha; };

  /**   \brief make an educated guess, what this function does...
  *
  */
  const Eigen::VectorXi &get_cw(void) const { return _cw; };

 protected:
  /**     \brief computes indices in weighted sparse grid space (recursive)
  *	    \param[in] maxBit position of the digit of the current multiindex
  *                  such that only digits with position greater or equal to
  *                  maxBit are considered in the recursion.
  *	 	 \param[in] k counter for the number of multiindices is passed
  *                  as pointer
  *	    \param[in] cpFun function/functor/whatever which decides whether a
  *                  multiindex is in the indexset or not.
  *
  *	    Additionally, the level of the sparse grid _q, the dimension of the
  *	    multiindices _dim and the memory _alpha for the indexset is
  *	    provided globally and referred to.
  *
  *	    The function Yw_alpha recursively computes the indexset of a
  *	    generalized sparse grid which is described by a level _q and a
  *	    function cpFun, i.e. it finds all multiindices alpha such that
  *	    cpFun(alpha)<= _q. To that end, we travel the set Xw_alpha and
  *       at a multi index only if the corresponding _cw is non zero.
  *	    The function needs to be called with input parameters
  *	    maxBit=0, k=0 defined as a reference, predefined
  *	    function cpFun which maps from \mathbb{N}^{_dim}\to \mathbb{R}.
  *	    Since cpFun(0)=0, the zero multiindex is always the first index
  *	    which is included in the sparse grid on default and serves as the
  *	    root node in the tree structure of the index set.
  *	    On level l of the tree, all indices alpha with
  *	    cpFun(alpha)<= _q and |alpha| = l-1 are computed.
  *	    Moreover, the tree is structured such that all indices alpha' which
  *	    are sons of a multiindex alpha fulfill alpha'>=alpha elementwise.
  *	    This structure immediately allows to exploit the downward closeness
  *	    of the indexset since we do not need a further recursion level
  *	    at an edge where the corresponding index does not belong to the
  *       indexset.
  *	    Moreover, the integer maxBit is necessary to avoid multiple
  *       repetition of multiindices. In every step maxBit is the position of
  *       the digit which is modified to obtain the current multiindex. In the
  *       subtree corresponding to this index only digits greater or equal to
  *       maxBit will be considered for modification. A simple example for the
  *       tree structure for _q=2, _dim=2 and cpFun(alpha) = alpha_1+alpha_2 is
  *       depicted below.
  *
  *								(0,0;mB=0;*k=0)
  *								/            \
  *						(1,0;mB=0;*k=1) (0,1;mB=1;*k=4)
  *						/		      \               \
  *			(2,0;mB=0;*k=2) (1,1;mB=1;*k=3)	(0,2;mB=1;*k=5)
  *
  */
  template <class Compare>
  void combiIndexSet(int maxBit, int *k, const Compare &cpFun,
                     Eigen::VectorXi &currInd) {
    int cw = 0;
    // successively increase all entries in the current index
    for (int i = maxBit; i < _dim; ++i) {
      ++currInd(i);
      if (cpFun(currInd) <= _q) {
        // if the index is in Xw_alpha, check if its weight is non zero
        if (cpFun(currInd + _myOnes) > _q)
          cw = combiWeights(0, 1, 1, currInd, cpFun);
        else
          cw = 0;
        // if its weight is non zero, add it to the index set
        if (cw) {
          if (_alpha.cols() <= *k) {
            _alpha.conservativeResize(_dim, _alpha.cols() + __MEMCHUNKSIZE__);
            _cw.conservativeResize(_cw.size() + __MEMCHUNKSIZE__);
          }
          _alpha.col(*k) = currInd;
          _cw(*k) = cw;
          ++(*k);
        }
        // check son indices only if father index is in Xw. this is possible
        // due to the downward closedness assumption
        combiIndexSet(i, k, cpFun, currInd);
      }
      --currInd(i);
    }
  }

  /**     \brief computes weights for the tensor product quadrature (recursive)
  *	    \param[in] maxBit position of the digit of the current multiindex
  *                  such that only digits with position greater to maxBit are
  *	 		         considered in the recursion.
  *       \param[in] cw is initialize by 1 due to alpha+0\in Xw if alpha in Xw
  *	    \param[in] lvl current level-1 in the recursion tree.
  *	    \param[in] ind current multi index alpha
  *	    \param[in] cpFun function/functor/whatever which decides whether a
  *                  multi index is in the indexset or not.
  *
  *	    Additionally, the level of the sparse grid _q, the dimension of the
  *	    multiindices _dim and the memory _alpha for the indexset is
  *	    provided globally.
  *
  *	    The function cw_alpha computes the coefficient of the indth
  *	    multi index in the combination technique formula corresponding to
  *	    the generalized sparse quadrature. This coefficient is given by
  *	    cw(alpha) = sum_{beta in {0,1}^_dim, alpha+beta in _alpha}
  *       (-1)^|beta|. Hence, in order to determine the value of the
  *       coefficient, we have to check for each beta in {0,1}^_dim if
  *       alpha+beta belongs to the index set _alpha. For beta=(0,...,0),
  *       alpha+beta belongs to _alpha and, hence, we
  *       initialize cw=1 and do not have to consider this multiindex.
  *	    The downward closeness of the index set can be exploited once more
  *	    which results in a similar tree structure as in the function
  *       Yw_alpha. The only difference of the tree is that once a digit of the
  *       multiindex alpha is modified, we do not have to consider further
  *       modification of this digit due to the condition beta in {0,1}^_dim.
  *       This is realized in the recursive algorithm consider only digits which
  *       are greater than maxBit in contrast to the "greater than and equal"
  *       condition in Yw_alpha.
  *	    In the recursive algorithm, the multiindices beta do not need to be
  *       computed explicitly, but the multiindex alpha is modified at each node
  *       of the tree. Of course, this modification has to be revoked when the
  *       recursion goes back to the parent node in the tree.
  *	    The value of lvl=|beta| is important since it determines whether the
  *	    coefficient is modified by +1 or -1.
  *	    A simple example of the algorithm for _q=2, _dim=2, ind = (1,0) and
  *	    cpFun(alpha) = alpha_1+alpha_2 is depicted below.
  *	    The further initial values are maxBit=mB=0, cw = 1, lvl=0.
  *
  *								(1,0;mB=0;lvl=0;cw=1)
  *								/                   \
  *					(2,0;mB=0;lvl=1;cw=0)   (1,1;mB=1;lvl=1;cw=-1)
  *							/
  *			[ (2,1;mB=1;lvl=2;cw=0) ] <-- is not contained in the
  *                                   indexset, no modification of cw
  *
  *	    The algorithm returns the final value of cw. In the example above
  *       this would be the correct value cw(1,0)=-1.
  */
  template <class Compare>
  int combiWeights(int maxBit, int cw, int lvl, Eigen::VectorXi &ind,
                   const Compare &cpFun) {
    // successively check all bits of the multiindex is contained
    for (int i = maxBit; i < _dim; ++i) {
      ++ind(i);
      if (cpFun(ind) <= _q) {
        if (lvl % 2)
          --cw;
        else
          ++cw;
        // again, exploit downward closedness and perform recursion only
        // if father is contained in Xw_alpha
        cw = combiWeights(i + 1, cw, lvl + 1, ind, cpFun);
      }
      --ind(i);
    }

    return cw;
  }

  // Member variables
  Eigen::MatrixXi _alpha;
  Eigen::VectorXi _cw;
  Eigen::VectorXi _myOnes;
  int _dim;
  int _q;
};
#endif
