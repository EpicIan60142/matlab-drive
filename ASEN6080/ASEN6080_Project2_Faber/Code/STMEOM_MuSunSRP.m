function dXPhi = STMEOM_MuSunSRP(t,XPhi,pConst,scConst)
% Function for the equations of motion for an STM including the sun as a
% 3rd body perturbation and SRP.
%   Inputs:
%       - t: Current integration time
%       - XPhi: State transition matrix and state vector (6^2 + 6)x1:
%               [X; Phi]
%       - pConst: Planetary constants structure as definted by
%                 getPlanetConst.m
%       - scConst: Spacecraft constants structure as defined by
%                  getSCConst.m
%   Outputs:
%       - dPhiX: Rate of change vector for the provided STM and state.
%                Note: dPhi is n^2x1, not nxn!
%                [dX; dPhi]
%
%   By: Ian Faber, 04/14/2025
%

x = XPhi(1);
y = XPhi(2);
z = XPhi(3);
xDot = XPhi(4);
yDot = XPhi(5);
zDot = XPhi(6);
C_R = XPhi(7);
Phi = XPhi(8:end);

n = 7;

Phi = reshape(Phi, n, n); % Phi is converted into a 49x1 vector...

X_a = [x; y; z; xDot; yDot; zDot; C_R];

r = sqrt(x^2 + y^2 + z^2);

A = DynamicsPartials_MuSunSRP(X_a, t, pConst, scConst);

dPhi_mat = A*Phi;
dPhi = reshape(dPhi_mat,n^2,1); % Need to convert back to a column vector...

X0 = [x;y;z;xDot;yDot;zDot;C_R];
dX = orbitEOM_MuSunSRP(t,X0,pConst,scConst);

dXPhi = [dX; dPhi];

end



