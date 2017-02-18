%
% Numerical Example 1 from
%
% A.-L. Haji-Ali, H. Harbrecht, M. Peters, and M. Siebenmorgen. 
% Novel results for the anisotropic sparse quadrature and their 
% impact on random diffusion problems. Preprint 2015-27, 
% Mathematisches Institut, Universität Basel, Switzerland, 2015.
%

clc;
clear all;
close all;
format long;

% decay parameter
r = 2;
% dimension
dim = 1000;
% maximum level q
maxLvl = 16;

% testfunction
fun = @(x) 1 ./ (0.6 + 0.2 * sum(kron(ones(1,size(x,2)),([1:dim].^-r)') .* x));

% reference solutions
ref2 = 1.7393632457035437;
ref3 = 1.7342253547471955;
ref4 = 1.7331866232415222;

% init cell array with univariate quadrature rules up to maxLvl
Quad = cell(maxLvl+1,1);
for i = 0:40
    [xi,w] = univariateGaussLegendre(i);
    Quad{i+1} = [xi;w];
end

% compute the weights for the sparse index set
kappa = 1./[1:dim].^-r + sqrt(1+1./[1:dim].^(-2*r));
w = log(kappa);

% evaluate the sparse quadrature
err = 0;
pts = 0;
for i = 0:maxLvl
tic
[Q,W,sort] = MXsparseQuadrature(i,dim,'TD', Quad, w);
quad = fun(Q) * W;
err(i+1) = abs(ref2-quad);
pts(i+1) = length(W);
display(sprintf('lvl: %3d numPts: %8d error: %8e', i, pts(i+1), err(i+1)));
toc
end

% plot error
figure(1);
loglog(pts,err,'k-s')
