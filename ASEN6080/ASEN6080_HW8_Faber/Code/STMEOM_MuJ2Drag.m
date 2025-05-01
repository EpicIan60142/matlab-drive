function dXPhi = STMEOM_MuJ2Drag(t,XPhi,pConst,scConst)
% Function for the equations of motion for an STM of the 2 body problem
% including J2 contributions
%   Inputs:
%       - t: Current integration time
%       - XPhi: State transition matrix and state vector (6^2 + 6)x1:
%               [X; Phi]
%       - mu: Gravitational parameter for the point mass of the central 
%             body of interest
%       - J2: Gravitational parameter for the oblateness of the central
%             body of interest
%       - Ri: Reference radius
%   
%   Outputs:
%       - dPhiX: Rate of change vector for the provided STM and state.
%                Note: dPhi is 6^2x1, not 6x6!
%                [dX; dPhi]
%
%   By: Ian Faber, 01/23/2025
%

x = XPhi(1);
y = XPhi(2);
z = XPhi(3);
xDot = XPhi(4);
yDot = XPhi(5);
zDot = XPhi(6);
Phi = XPhi(7:end);

Phi = reshape(Phi, 6, 6); % Phi is converted into a 49x1 vector...

X_a = [x; y; z; xDot; yDot; zDot];

r = sqrt(x^2 + y^2 + z^2);

A = DynamicsPartials_MuJ2Drag(X_a, pConst, scConst);

dPhi_mat = A*Phi;
dPhi = reshape(dPhi_mat,36,1); % Need to convert back to a column vector...

X0 = [x;y;z;xDot;yDot;zDot];
dX = orbitEOM_MuJ2Drag(t,X0,pConst,scConst);

dXPhi = [dX; dPhi];

end



