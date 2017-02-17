#include <string>
#include <iostream>

#include <Eigen/Dense>

#include "mex.h"
#include "matrix.h"

#include "SparseIndexSet.hpp"
#include "TDindexSet.hpp"
#include "HCindexSet.hpp"
#include "SparseQuadrature.hpp"
#include "Cell2univariateQuadrature.hpp"

// fancy compare functor that performs the callback to Matlab
struct MatlabCpFun {
  const mxArray *_prhs[2];

  mxArray *_MatAlpha;
  int _dim;
  MatlabCpFun(const mxArray *prhs, int dim) {
    _MatAlpha = mxCreateDoubleMatrix(dim, 1, mxREAL);
    _prhs[0] = prhs;
    _prhs[1] = _MatAlpha;
    _dim = dim;
  };
  double operator()(const Eigen::VectorXi &alpha) const {
    mxArray *plhs;
    Eigen::Map<Eigen::VectorXd>(mxGetPr(_MatAlpha), _dim) =
        alpha.cast<double>();
    mexCallMATLAB(1, &plhs, 2, const_cast<mxArray **>(_prhs),
                  (const char *)"feval");
    if (mxIsScalar(plhs))
      return *(mxGetPr(plhs));
    else
      mexErrMsgIdAndTxt(
          "MATLAB:MXsparseQuadraturecpp",
          "MXsparseQuadrature requires adequate function handle. Check help.");
    return 0.;
  };
};

/** nlhs Number of expected output mxArrays
*   plhs Array of pointers to the expected output mxArrays
*   nrhs Number of input mxArrays
*   prhs Array of pointers to the input mxArrays.
*        Do not modify any prhs values in your MEX file.
*        Changing the data in these read-only mxArrays can
*        produce undesired side effects.
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int q = 0;
  int dim = 0;
  int nPts = 0;
  int myCASE = 0;
  TDindexSet TDind;
  HCindexSet HCind;
  SparseIndexSet GENind;
  SparseQuadrature Q;
  Cell2univariateQuadrature CellUniQ;

  if (nrhs != 5) {
    mexErrMsgIdAndTxt(
        "MATLAB:MXsparseQuadraturecpp:nargin",
        "MXsparseQuadrature requires five input arguments. Check help.");
  } else if (nlhs < 2) {
    mexErrMsgIdAndTxt("MATLAB:MXsparseQuadraturecpp:nargin",
                      "MXsparseQuadrature requires at least two output "
                      "arguments. Check help.");
  } else if (!mxIsScalar(prhs[0]) || !mxIsScalar(prhs[1])) {
    mexErrMsgIdAndTxt("MATLAB:MXsparseQuadraturecpp",
                      "First two arguments have to be scalar. Check help.");
  } else if (!mxIsChar(prhs[2])) {
    mexErrMsgIdAndTxt("MATLAB:MXsparseQuadraturecpp",
                      "Third argument has to be a string. Check help.");
  } else if (!mxIsCell(prhs[3]))
    mexErrMsgIdAndTxt("MATLAB:MXsparseQuadraturecpp",
                      "Fourth argument has to be a Cell array. Check help.");

  q = std::round(*(mxGetPr(prhs[0])));
  dim = std::round(*(mxGetPr(prhs[1])));
  std::string type(mxArrayToString(prhs[2]));

  CellUniQ.initQuadrature(prhs[3]);
  
  // convert type to uppercase
  for (auto i = type.begin(); i != type.end(); ++i)
    *i = std::toupper(*i);

  if (type == "TD") {
    if (mxIsClass(prhs[4], "function_handle"))
      mexErrMsgIdAndTxt(
          "MATLAB:MXsparseQuadraturecpp",
          "Last argument has to be array if TD is used. Check help.");
    TDind.computeIndexSet(q, Eigen::Map<Eigen::VectorXd>(mxGetPr(prhs[4]), dim));
    Q = SparseQuadrature(TDind, CellUniQ);
    Q.purgeSparseQuadrature();
    plhs[2] = mxCreateDoubleMatrix(dim, 1, mxREAL);
    Eigen::Map<Eigen::VectorXd> sort(mxGetPr(plhs[2]), dim);
    const Eigen::VectorXi &mySort = TDind.get_sortW();
    sort = mySort.cast<double>().array() + 1;
  } else if (type == "HC") {
    if (mxIsClass(prhs[4], "function_handle"))
      mexErrMsgIdAndTxt(
          "MATLAB:MXsparseQuadraturecpp",
          "Last argument has to be array if HC is used. Check help.");

    HCind.computeIndexSet(q, Eigen::Map<Eigen::VectorXd>(mxGetPr(prhs[4]), dim));
    Q = SparseQuadrature(HCind, CellUniQ);
    Q.purgeSparseQuadrature();
  } else if (type == "GEN") {
    if (!mxIsClass(prhs[4], "function_handle"))
      mexErrMsgIdAndTxt("MATLAB:MXsparseQuadraturecpp",
                        "Last argument has to be function handle if Gen is "
                        "used. Check help.");

    GENind.computeIndexSet(q, dim, MatlabCpFun(prhs[4], dim));
    Q = SparseQuadrature(GENind, CellUniQ);
    Q.purgeSparseQuadrature();
  } else {
    mexErrMsgIdAndTxt("MATLAB:MXsparseQuadraturecpp",
                        "type has an invalid value. Check help.");
  }

  const Eigen::MatrixXd &myqPoints = Q.get_qPoints();
  const Eigen::VectorXd &myqWeights = Q.get_qWeights();

  nPts = myqWeights.size();

  plhs[0] = mxCreateDoubleMatrix(dim, nPts, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(nPts, 1, mxREAL);

  if (nlhs == 3 && type != "TD") plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);

  // set quadrature Points and Weights
  Eigen::Map<Eigen::MatrixXd> qPoints(mxGetPr(plhs[0]), dim, nPts);
  Eigen::Map<Eigen::VectorXd> qWeights(mxGetPr(plhs[1]), nPts);

  qPoints = myqPoints;
  qWeights = myqWeights;

  return;
}
