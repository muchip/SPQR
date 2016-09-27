#include "SparseQuadrature.hpp"

/**   \brief non-void constructor
*     \param[in] spInd is an object of the class sparseIndexSet
*     \param[in] Q is an object of the class univariateQuadrature
*
*/
SparseQuadrature::SparseQuadrature(const SparseIndexSet &spInd,
                                   UnivariateQuadrature &Q) {
  SparseQuadrature::computeSparseQuadrature(spInd, Q);
}

/**   \brief initializes the sparse quadrature based on the sparse index set
*            and the tensor product quadrature
*     \param[in] spInd is an object of the class sparseIndexSet
*     \param[in] Q is an object of the class univariateQuadrature
*
*/
void SparseQuadrature::computeSparseQuadrature(const SparseIndexSet &spInd,
                                               UnivariateQuadrature &Q) {
  int nPts = 0;
  int M = 0;
  // get references to the index set alpha and the corresponding weights cw
  const Eigen::MatrixXi &alpha = spInd.get_alpha();
  const Eigen::VectorXi &cw = spInd.get_cw();
  M = (int)alpha.rows();
  // resize the arrays for quadrature points and weights
  _qPoints.resize(M, __MEMCHUNKSIZE__);
  _qWeights.resize(__MEMCHUNKSIZE__);
  // initialize all tensor product quadrature points and weights
  for (int j = 0; j < alpha.cols(); ++j) {
    TensorProductQuadrature locTPQ(alpha.col(j), Q);
    // resize if necessary
    if (_qPoints.cols() <= nPts + locTPQ.get_nPts()) {
      if (locTPQ.get_nPts() < __MEMCHUNKSIZE__) {
        _qPoints.conservativeResize(M, _qPoints.cols() + __MEMCHUNKSIZE__);
        _qWeights.conservativeResize(_qWeights.size() + __MEMCHUNKSIZE__);
      } else {
        _qPoints.conservativeResize(M, _qPoints.cols() + locTPQ.get_nPts());
        _qWeights.conservativeResize(_qWeights.size() + locTPQ.get_nPts());
      }
    }
    _qPoints.block(0, nPts, M, locTPQ.get_nPts()) = locTPQ.get_points();
    _qWeights.segment(nPts, locTPQ.get_nPts()) = cw(j) * locTPQ.get_weights();
    // add up total number of quadrature points
    nPts += locTPQ.get_nPts();
  }
  // crop arrays to the correct size
  _qPoints.conservativeResize(M, nPts);
  _qWeights.conservativeResize(nPts);
}

/**   \brief removes duplicate points from the quadrature and combines their
*            weights, i.e. reduce combination technique to the corresponding
*            generalized sparse grid quadrature.
*            application may result in a major speedup if the evaluation of
*            each quadrature point is expensive, since it avoids double
*            evaluations of a given point
*
*/
void SparseQuadrature::purgeSparseQuadrature(void) {
  // get variables to store sorted arrays
  Eigen::VectorXi sortVec;
  Eigen::VectorXi uniquePoints;
  double dtmp = 0;
  int itmp = 0;
  int jtmp = 0;
  int nPts = 0;
  int M = 0;
  int currind = 0;
  // get array sizes
  nPts = (int)_qPoints.cols();
  M = (int)_qPoints.rows();
  // initialize the index vector 0:1:nPts-1
  sortVec = Eigen::ArrayXi::LinSpaced(nPts, 0, nPts - 1);
  uniquePoints.resize(nPts);
  uniquePoints.setZero();
  // use the C++ quick sort to sort the array _qPoints lexicographically
  // this is reflected by the permutations in sortVec cost is in avarage
  // nPts*log(nPts)*M, where we use the pessimistic estimate M for the
  // lexicographical comparison
  std::sort(&(sortVec(0)), &(sortVec(0)) + nPts, lexiCompareInd(_qPoints));
  // now, we update the weights of the sorted quadrature points according to
  // their multiplicity the indices of the unique quadrature points together
  // with their respective weights are found in uniquePoints
  currind = 0;
  uniquePoints(currind) = sortVec(0);
  // this loop checks if the next point is equal to the previous one. If true,
  // the weights are added else the point index is added to uniquePoints
  for (int i = 1; i < nPts; ++i) {
    if (!isEqual(_qPoints.col(uniquePoints(currind)),
                 _qPoints.col(sortVec(i)))) {
      // check if summed up quadrature weight is significant else dump the
      // current point
      if (std::abs(_qWeights(uniquePoints(currind))) > __PRECISION__) ++currind;
      uniquePoints(currind) = sortVec(i);
    } else
      _qWeights(uniquePoints(currind)) += _qWeights(sortVec(i));
  }
  // crop uniquePoints and sortVec to the correct length
  // now use sortVec to get access to the inverse indices in uniquePoints
  uniquePoints.conservativeResize(currind + 1);
  sortVec.setZero();
  for (int i = 0; i <= currind; ++i) sortVec(uniquePoints(i)) = i;
  // now, move everything in place
  for (int i = 0; i < uniquePoints.size(); ++i)
    if (uniquePoints(i) != i) {
      itmp = uniquePoints(i);
      jtmp = sortVec(i);
      _qPoints.col(i).swap(_qPoints.col(itmp));
      dtmp = _qWeights(i);
      _qWeights(i) = _qWeights(itmp);
      _qWeights(itmp) = dtmp;
      uniquePoints(i) = i;
      sortVec(i) = i;
      uniquePoints(jtmp) = itmp;
      sortVec(itmp) = jtmp;
    }
  _qPoints.conservativeResize(M, currind + 1);
  _qWeights.conservativeResize(currind + 1);
}

/**   \brief make an educated guess, what this function does...
*
*/
const Eigen::MatrixXd &SparseQuadrature::get_qPoints(void) const {
  return _qPoints;
}

/**   \brief make an educated guess, what this function does...
*
*/
const Eigen::VectorXd &SparseQuadrature::get_qWeights(void) const {
  return _qWeights;
}
