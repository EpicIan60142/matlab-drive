function dX = calculateAttitude(X, I, u)
% Calculates the nano-satellite body attitude relative to inertial space
%
%   Inputs:
%       - X: State vector at a given point in time
%               [sig_1; sig_2; sig_3; w_1; w_2; w_3]
%       - I: Body fixed inertia matrix
%               [diag(I_11, I_22, I_33)]
%       - u: Control input vector
%               [ u_1; u_2; u_3]
%   Outputs:
%       - dX: Rate of change vector based on the current state
%               [sigDot; wDot]
%

sig = X(1:3);
w = X(4:6);

sigSqr = dot(sig,sig);
sigDot = 0.25*((1-sigSqr)*eye(3) + 2*tilde(sig) + 2*(sig*sig'))*w;

wDot = (I^-1)*(-tilde(w)*I*w + u);

dX = [sigDot; wDot];

end