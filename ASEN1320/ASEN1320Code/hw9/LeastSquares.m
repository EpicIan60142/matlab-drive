function [m,b] = LeastSquares(x,y)

N = length(x);

A = 0; B = 0; C = 0; D = 0;

for j = 1:N
    A = A + x(j);
    B = B + y(j);
    C = C + x(j)*y(j);
    D = D + x(j)^2;
end

m = (A*B - N*C)/(A^2 - N*D);
b = (A*C - B*D)/(A^2 - N*D);
 
end
