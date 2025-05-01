function out = Householder(mat)
% Function that implements the Householder Algorithm on a given matrix
%   Inputs:
%       - mat: Generic (m+n)x(n+1) matrix of the form
%              [Rbar, bbar; H, y]
%   Outputs:
%       - out: (m+n)x(n+1) Upper triangularized matrix of the form
%              [R, b; 0, e]
%
%   By: Ian Faber, 03/07/2025
%

%% Make copy of input matrix
A = mat;

%% Find sizes
n = size(A,2) - 1;
m = size(A,1) - n;

%% Run algorithm
out = zeros(size(A));
for k = 1:n
        % Pull out main diagonal element and the elements of its column at
        % and below the diagonal
    A_kk = A(k,k);
    A_ik = A(k:m+n,k);

        % Compute sigma - the norm of the original vector, and u_k - 
        % the first element of the reflection vector
    sigma = sign(A_kk)*sqrt(sum(A_ik.^2));
    u_k = A_kk + sigma;

        % Apply transformation to main diagonal and assign the elements
        % of the reflection vector
    A(k,k) = -sigma;
    A_ik = A(k+1:m+n,k);
    u_i = [u_k; A_ik];

        % Calculate beta
    beta = 1/(sigma*u_k);

        % Apply transformation to every other column
    for j = k+1:n+1
        A_ij = A(k:m+n,j);
        
        gamma = beta*sum(u_i.*A_ij);

        A(k:m+n,j) = A_ij - gamma*u_i;
    end

        % Zero out everything under the main diagonal
    A(k+1:m+n,k) = zeros(size(A(k+1:m+n,k)));

end

out = A;

end