clc;
clear all;
close all;
format long;

fun1 = @(x) prod((1.+1./size(x,1))*(0.5*x+0.5).^(1./size(x,1)),1);
fun2 = @(x) prod(0.5*pi*sin(pi*(0.5*x+0.5)),1);
fun3 = @(x) prod(exp(0.5*x+0.5)/(exp(1)-1));
maxLvl = 4;
dim = 3;


% init cell array with univariate quadrature rules up to maxLvl
Quad = cell(maxLvl+1,1);
for i = 0:maxLvl
    [xi,w] = univariateClenshawCurtis(i);
    %[xi,w] = univariateGaussLegendre(i);
    %[xi,w] = univariateTrapezoidalRule(i);
    Quad{i+1} = [xi;w];
end

% test sparse quadrature
for i = 0:maxLvl
tic
[Q,W,sort] = MXsparseQuadrature(i,dim,'Gen', Quad, @(x) sum(x));
%[Q,W,sort] = MXsparseQuadrature(i,dim,'TD', Quad, ones(1,dim));
%[Q,W,sort] = MXsparseQuadrature(i,dim,'HC', Quad, ones(1,dim));
quad = fun1(Q) * W;
err = abs(1-quad);
display(sprintf('lvl: %3d numPts: %8d error: %8e', i, length(W), err));
toc
end
