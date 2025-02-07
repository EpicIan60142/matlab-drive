function val = compSimp(x, f)
% compSimp - Composite simpson's rule function
%       Implements composite simpson's rule formula
%       Inputs:
%           x: vector of independent variables
%           f: function handle of function to integrate
%       Outputs:
%           val: resultant value of integration
%       
%       Author: Ian Faber
%       Collaborators: Maggie Wussow
%       Date: 1/28/2023

% n = 2N, n + 1 points in x

val = [];

a = x(1);
b = x(end);

n = length(x) - 1;

h = (b-a)/n;

for k = 1:(n/2)
    f_1 = f(x(2*k - 1));
    f_2 = f(x(2*k));
    f_3 = f(x(2*k + 1));
    val = [val; f_1 + 4*f_2 + f_3];
end

val = (h/3)*sum(val);

end