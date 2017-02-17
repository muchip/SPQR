function [xi,w] = univariateTrapezoidalRule(l)
% abscissas and weights for the nested open trapezoidal rule, I=[-1,1]
% the weight is normalized (probability measure)
% l is the quadrature level

wRef = [0.5 0.5];
xiRef = [1/3 2/3];

tau = @(x,l,i) 2^(1-l)*(x+i);

n = 2^l;

if n == 1
    % midpointrule for the first Level
    xi = 0;
    w = 2;
else
    xi = zeros(1,n);
    for j=1:2^l
        xi(2*(j-1)+1:2*(j-1)+2) = tau(xiRef,l+1,j-1);
    end
    xi = 2*xi-1;
    w = ones(size(xi))*0.5/2^(l);
end
end
