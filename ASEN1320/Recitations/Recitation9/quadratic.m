function [x1,x2] = quadratic(a, b, c)
% Computes the roots of the quadratic equation
% Inputs: a, b, c - coefficients
% Outputs: x1, x2 - roots

x1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
x2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);

%Don't need a return statement, just assign values to outputs!

end

