%   MXsparseQuadrature.m Help file for sparse Quadrature MEX-file.
%
%   [Q,W,sort] = MXsparseQuadrature(q,dim,type,Quad,w/CpFun);
%
%   INPUT:
%   q     maximum level in the sparse grid quadrature
%   dim   dimension of the sparse grid quadrature
%   type  either 'HC' for hyperbolic cross, 'TD' for total degree
%         or 'Gen' for a general sparse grid provided by the criterion
%         in the fourth argument CpFun as function handle
%   Quad  contains a cell array with univariate quadrature rules
%         Quad{i}=[xi;w], where xi and w are row vectors
%         in the weighted case, Quad should contain the rules for levels
%         0..q/max(w)
%   w/CpFun either weight vector for 'HC' and 'TD' or function handle
%           to determine for a multiindex CpFun(alpha) <= q
%
%   OUTPUT:
%   Q     contains the quadrature points
%   W     contains the quadrature weights
%   sort  if 'TD' is used, contains the order of the dimensions in concordance
%         with the increasingly sorted vector w. is not referred to for 'HC' and
%         'Gen'
