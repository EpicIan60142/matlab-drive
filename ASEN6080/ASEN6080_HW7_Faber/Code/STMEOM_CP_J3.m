function dXTheta = STMEOM_CP_J3(t, XTheta, pConst)
% Function for the equations of motion for an STM of the 2 body problem
% including J2 contributions with J3 as a consider parameter
%   Inputs:
%       - t: Current integration time
%       - XTheta: State transition matrix and state vector (6 + 6)x1:
%                 [X; Theta]
%       - pConst: Planetary constants structure as defined in
%                 getPlanetConst.m
%   
%   Outputs:
%       - dThetaX: Rate of change vector for the provided STM and state.
%                  Note: dTheta is 6x1!
%                  [dX; dTheta]
%
%   By: Ian Faber, 04/06/2025
%

x = XTheta(1);
y = XTheta(2);
z = XTheta(3);
xDot = XTheta(4);
yDot = XTheta(5);
zDot = XTheta(6);
Theta = XTheta(7:12);

% Theta = reshape(Theta, 6, 1);

X_a = [x; y; z; xDot; yDot; zDot];

r = sqrt(x^2 + y^2 + z^2);

A = DynamicsPartials_MuJ2J3(X_a, pConst.mu, pConst.J2, pConst.J3, pConst.Ri);

B = DynamicsPartials_CP_J3(X_a, pConst);

dTheta = A*Theta + B;
% dTheta = reshape(dTheta_mat,6,1);

X0 = [x; y; z; xDot; yDot; zDot];
dX = orbitEOM_MuJ2(t,X0,pConst.mu,pConst.J2,pConst.Ri);

dXTheta = [dX; dTheta];

end



