function dXPhi = STMEOM_MuJ2Drag(t,XPhi,pConst,scConst)
% Function for the equations of motion for an STM of the 2 body problem
% including Mu, J2, and drag acceleration contributions
%   Inputs:
%       - t: Current integration time
%       - XPhi: State transition matrix and state vector (18^2 + 18)x1:
%               [X; Phi]
%       - pConst: Planetary constants struct as defined by getPlanetConst.m
%       - scConst: Spacecraft constants struct as defined by getSCConst.m
%   
%   Outputs:
%       - dPhiX: Rate of change vector for the provided STM and state.
%                Note: dPhi is 6^2x1, not 6x6!
%                [dX; dPhi]
%
%   By: Ian Faber, 02/11/2025
%

x = XPhi(1);
y = XPhi(2);
z = XPhi(3);
xDot = XPhi(4);
yDot = XPhi(5);
zDot = XPhi(6);
mu = XPhi(7);
J2 = XPhi(8);
Cd = XPhi(9);
Phi = XPhi(19:end);

Phi = reshape(Phi, 18, 18); % Phi is converted into a column vector...

X = XPhi(1:18);

A = DynamicsPartials_MuJ2Drag(X, pConst, scConst);

dPhi_mat = A*Phi;
dPhi = reshape(dPhi_mat,18^2,1); % Need to convert back to a column vector...

dX = orbitEOM_MuJ2Drag(t,X,pConst,scConst);

dXPhi = [dX; dPhi];

end



