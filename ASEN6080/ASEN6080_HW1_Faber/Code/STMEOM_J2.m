function dPhiX = STMEOM_J2(t,XPhi,mu,Ri)
% Function for the equations of motion for an STM of the 2 body problem
% including J2 contributions
%   Inputs:
%       - t: Current integration time
%       - XPhi: State transition matrix and state vector (7^2 + 7)x1:
%               [X; Phi]
%       - mu: Gravitational parameter for the central body of interest
%       - Ri: Reference radius
%   
%   Outputs:
%       - dPhiX: Rate of change vector for the provided STM and state.
%                Note: dPhi is 7^2x1, not 7x7!
%                [dPhi; dX]
%
%   By: Ian Faber, 01/23/2025
%

x = XPhi(1);
y = XPhi(2);
z = XPhi(3);
xDot = XPhi(4);
yDot = XPhi(5);
zDot = XPhi(6);
J2 = XPhi(7);
Phi = XPhi(8:56);

Phi = reshape(Phi, 7, 7); % Phi is converted into a 49x1 vector...

X_a = [x; y; z; xDot; yDot; zDot; mu; J2];

r = sqrt(x^2 + y^2 + z^2);

A = DynamicsPartials_MuJ2(X_a, Ri);

dPhi_mat = A*Phi;
dPhi = reshape(dPhi_mat,49,1); % Need to convert back to a column vector...

X0 = [x;y;z;xDot;yDot;zDot;J2];
dX = orbitEOM_MuJ2(t,X0,mu,Ri);

dPhiX = [dX; dPhi];

end



