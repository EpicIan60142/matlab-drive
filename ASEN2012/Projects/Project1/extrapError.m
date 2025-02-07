function error = extrapError(W, A, t)
% Function that calculates the extrapolated error of a regression line at
% some new time index
%   Inputs: Weight matrix W, Least Squares matrix A, time index t
%   Outputs: Extrapolated regression error

Q = (A'*W*A)^-1;
error = sqrt([t 1]*Q*[t; 1]);

end

