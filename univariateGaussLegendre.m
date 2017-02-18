function [xi,w] = univariateGaussLegendre (l)
% abscissas and weights for the Gauss Legendre formulas, I=[-1,1]
% the weight is normalized (probability measure)
% l is the level of the quadrature

n = ceil(0.5*l) + 1;

if n == 1
    % midpointrule for the first Level
    xi = 0;
    w = 1;
else
   subd = [1:n-1] ./ sqrt(4*[1:n-1].^2-1);
   A = zeros(n,n);
   A = diag([subd'],-1) + diag([subd'],1);
   [V,D] = eig(A);
   xi = diag(D)';
   w = V(1,:).^2;
   if mod(n, 2)
       xi(floor(n/2)+1) = 0;
end
end
