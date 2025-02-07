function x_hat = linReg(A, b)
% Function that carries out a least squares linear regression and outputs
% the slope and intercept of the regression
%   Inputs: Least squares matrix A, solution vector b
%   Outputs: Vector x_hat composed of the regression slope and intercept
    
x_hat = (A'*A)^-1*A'*b;

end

