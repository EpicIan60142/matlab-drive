function dX = rigidBodyEOM_MRP(X, L, I)
% Function that defines the EOM for a rigid body in terms of MRPs
%   Inputs:
%       - X: State vector at a given point in time
%               [sig1; sig2; sig3; w1; w2; w3]
%       - L: Control torque at a given point in time
%               [L1; L2; L3]
%       - I: Inertia matrix of the rigid body
%
%   Outputs:
%       - dX: Rate of change vector based on the current state
%               [dSig1; dSig2; dSig3; dw1; dw2; dw3]
%
    sig = X(1:3);
    w = X(4:6);

    sigSquared = norm(sig)^2;

    dSigma = 0.25*((1-sigSquared)*eye(3) + 2*tilde(sig) + 2*(sig*sig'))*w;
    dw = (I^-1)*(-tilde(w)*I*w + L);

    dX = [dSigma; dw];
end