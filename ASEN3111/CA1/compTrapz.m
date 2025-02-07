function val = compTrapz(x, f)
% compTrapz - Composite trapezoidal rule function
%       Implements the composite trapezoidal rule formula
%       Inputs: 
%           x: vector of independent variables 
%           f: function handle of function to integrate 
%       Outputs: 
%           val: resultant value of integration
%
%       Author: Ian Faber
%       Date: 1/21/2023

val = sum(diff(x).*movsum(f(x), 2, 'EndPoints', 'discard')/2);

end