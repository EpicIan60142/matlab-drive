% syntax: function [outputs] = name(inputs) ... end
% note that comments on the function purpose, inputs, and outputs can be
% helpful

function [x1,x2] = quadratic(a,b,c)
% Returns the roots of the quadratic equation
% Inputs:   a,b,c - first, second, third coefficients respectively
% Outputs:  x1,x2 - first and second roots respectively

x1 = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
x2 = (-b - sqrt(b^2 - 4*a*c)) / (2*a);

end