function [xi,w] = univariateClenshawCurtis (l)
% abscissas and weights for the nested Clenshaw-Curtis formulas, I=[-1,1]
% the weight is normalized (probability measure)
% abscissas are the extreme points of the Chebyshev polynomials
% l is the level of the quadrature

if l == 0
    n = 1;
else
    n = 2^l+1;
end

if n == 1
    % midpointrule for the first Level
    xi = 0;
    w = 1;
else
    ivec = 1:n;
    % evaluation points
    theta = (pi*(ivec-1)) / (n-1);
    % abscissa in the interval [-1,1]
    xi = cos(theta);
    % set the zero explicitly
    xi(2^(l-1)+1) = 0;
    % compute inner weights
    kvec = 1:n-2;
    % 2*cos(pi*(k-1))+2 = 4 for k-1 even,
    lambda = 4 ./ (-(kvec-1).^3-3*(kvec-1).^2+(kvec-1)+3);   
    % 2*cos(pi*(k-1))+2 = 0 for k-1 odd,
    lambda(2:2:end) = 0;   
    % discrete sine transform I
    lambda = dst(lambda);         
    w = 2*sin(theta(kvec+1)) .* lambda ./ ((n-1)*(1-xi(kvec+1).^2));
    % compute first and last weight w_1 and w_n
    gamma = (cos(pi*(ivec-1))+1) ./ (1-(ivec-1).^2);
    % prevent 0/0 situation for s = 1
    gamma(2) = 0;
    w1 = (gamma(1)+2*sum(gamma(2:n-1))+gamma(n)) / (2*n-2);
    w = 0.5 * [w1,w,w1];
end
end
