function [x1, x2] = quadratic(a, b, c)
% computes the roots of the quadratic equation
% inputs: a,b,c - first, second, third coefficents
% outputs: x1,x2 - first and second roots

x1 = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
x2 = (-b - sqrt(b^2 - 4*a*c)) / (2*a);

end

