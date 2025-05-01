function dXPhi = STMEOM_MuSunSRP_DMC(t,XPhi,B,u,pConst,scConst)
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
wX = XPhi(8);
wY = XPhi(9);
wZ = XPhi(10);
Phi = XPhi(11:end);

n = 10;

Phi = reshape(Phi, n, n); % Phi is converted into a 49x1 vector...

X_a = [x; y; z; xDot; yDot; zDot; C_R; wX; wY; wZ];

r = sqrt(x^2 + y^2 + z^2);

A = DynamicsPartials_MuSunSRP_DMC(X_a, t, B, pConst, scConst);

dPhi_mat = A*Phi;
dPhi = reshape(dPhi_mat,n^2,1); % Need to convert back to a column vector...

X0 = [x;y;z;xDot;yDot;zDot;C_R;wX;wY;wZ];
dX = orbitEOM_MuSunSRP_DMC(t,X0,B,u,pConst,scConst);

% C = [zeros(6,3); eye(3); zeros(1,3)];
% dX = dX + 1*C*u;

dXPhi = [dX; dPhi];

end



