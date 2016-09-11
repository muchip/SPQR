#ifndef __CELL2UNIVARIATEQUADRATURE__CLASS__
#define __CELL2UNIVARIATEQUADRATURE__CLASS__

#include <iostream>
#include <iomanip>
#include "univariateQuadrature.hpp"

#include <Eigen/Dense>

// include Matlab headers
#include "mex.h"
#include "matrix.h"

/** \brief extends univariate quadrature and wraps univariate quadrature
*          formulas provided by matlab into a univariateQuadrature
*/
class cell2univariateQuadrature : public univariateQuadrature {
 public:
  void initQuadrature(int maxLvl){
      // dummy routine to satisfy compiler, since virtual method has
      // to be instantiated
  };

  /** \brief transforms cell array into univariate quadrature some basic error
  *          checking is done. If an error occurs, Matlab takes care of it.
  */
  void initQuadrature(const mxArray *cell) {
    int n = 0;
    int m = 0;
    int ldCell = 0;
    const mxArray *pCell = NULL;
    const mwSize *dims = NULL;
    mwSize numDims = 0;

    numDims = mxGetNumberOfDimensions(cell);
    dims = mxGetDimensions(cell);
    ldCell = dims[0] < dims[1];

    if (numDims != 2 || dims[1 - ldCell] != 1)
      mexErrMsgIdAndTxt(
          "MATLAB:cell2univariateQuadrature",
          "cell array has to be of size (1,maxLvl+1) or (maxLvl+1,1).");

    // get memory for the quadrature vector
    _maxLvl = dims[ldCell] - 1;
    _Q.resize(_maxLvl + 1);
    for (mwSize i = 0; i <= (mwSize)_maxLvl; ++i) {
      // get pointer to cell entry and determine format
      pCell = mxGetCell(cell, i);
      m = mxGetM(pCell);
      n = mxGetN(pCell);
      if (m != 2)
        mexErrMsgIdAndTxt(
            "MATLAB:cell2univariateQuadrature",
            "each entry in the cell array has to be of the form [xi; "
            "w] where xi and w are row vectors.");
      // write points and weights to quadrature rule
      _Q[i].xi.resize(n);
      _Q[i].w.resize(n);
      double *Array = mxGetPr(pCell);
      for (int j = 0; j < n; ++j) {
        _Q[i].xi(j) = Array[2 * j];
        _Q[i].w(j) = Array[2 * j + 1];
      }
    }
  };

  /**    \brief checks if quadrature rules up to degree maxdeg exist
  *             if not, returns an error
  *      \param[in] maxLvl  maximum levl of the quadrature which is required
  */
  void resizeQuadrature(int maxLvl) {
    if (_maxLvl < maxLvl)
      mexErrMsgIdAndTxt("MATLAB:cell2univariateQuadrature",
                        "Too few univariate quadrature rules provided.");
  };

  /**    \brief testing routine for the Cell to univariate Quadrature.
  *      \param[in] maxLvl quadratures from lvl=0..maxLvl are printed
  *                 if maxLvl <= _maxLvl
  */
  void testQuadrature(int maxLvl) {
    maxLvl = _maxLvl < maxLvl ? _maxLvl : maxLvl;
    double sumw = 0;
    for (int i = 0; i <= maxLvl; ++i) {
      sumw = 0;
      for (int j = 0; j < _Q[i].xi.size(); ++j) {
        std::cout << std::setprecision(6) << "xi = " << std::setw(9)
                  << _Q[i].xi(j) << "\t"
                  << "w= " << std::setw(9) << _Q[i].w(j) << std::endl;
        sumw += _Q[i].w(j);
      }
      std::cout << "\nsum of weights: " << sumw << std::endl << std::endl;
      (void)maxLvl;
    };
  };
};
#endif
