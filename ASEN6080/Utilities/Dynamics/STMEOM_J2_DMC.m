function dXPhi = STMEOM_J2_DMC(t,XPhi,B,u,mu,J2,Ri)
% Function for the equations of motion for an STM of the 2 body problem
% including J2 contributions, including dynamic model compensation
%   Inputs:
%       - t: Current integration time
%       - XPhi: State transition matrix and state vector (6^2 + 6)x1:
%               [X; Phi]
%       - B: DMC time constant matrix organized as a 3x3 matrix:
%            B = diag([tau_x^-1, tau_y^-1, tau_z^-1])
%       - u: Gaussian white process noise organized as a 3x1 vector:
%            u = [uX; uY; uZ]
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
wX = XPhi(7);
wY = XPhi(8);
wZ = XPhi(9);
Phi = XPhi(10:end);

Phi = reshape(Phi, 9, 9); % Phi is converted into a 81x1 vector...

X_a = [x; y; z; xDot; yDot; zDot; wX; wY; wZ];

r = sqrt(x^2 + y^2 + z^2);

Aprime = DynamicsPartials_MuJ2_DMC(X_a, B, mu, J2, Ri);

dPhi_mat = Aprime*Phi;
dPhi = reshape(dPhi_mat,81,1); % Need to convert back to a column vector...

X0 = [x;y;z; xDot;yDot;zDot; wX;wY;wZ];
dX = orbitEOM_MuJ2_DMC(t,X0,B,mu,J2,Ri);

C = [zeros(6,3); eye(3)];
dX = dX + 0*C*u;

dXPhi = [dX; dPhi];

end



