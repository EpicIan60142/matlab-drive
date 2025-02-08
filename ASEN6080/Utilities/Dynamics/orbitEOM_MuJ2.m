function dX = orbitEOM_MuJ2(t,X,mu,J2,Ri)
% Function for the equations of motion for a body in orbit around a central
% body including mu and J2 in the 2 body problem
%   Inputs:
%       - t: Current integration time
%       - X: State vector arranged as follows:
%            [X; Y; Z; Xdot; Ydot; Zdot]
%       - mu: Gravitational parameter for the point mass of the central 
%             body of interest
%       - J2: Gravitational parameter for the oblateness of the central
%             body of interest
%       - Ri: Reference radius
%   
%   Outputs:
%       - dX: Rate of change vector for the provided state
%             [Xdot; Ydot; Zdot; Xddot; Yddot; Zddot]
%
%   By: Ian Faber, 01/23/2025
%

x = X(1);
y = X(2);
z = X(3);
xDot = X(4);
yDot = X(5);
zDot = X(6);

r = sqrt(x^2 + y^2 + z^2);

xDDot = -(mu*x/(r^3))*(1 + (((Ri/r)^2)*J2*((15/2)*(z/r)^2-(3/2))));
yDDot = -(mu*y/(r^3))*(1 + (((Ri/r)^2)*J2*((15/2)*(z/r)^2-(3/2))));
zDDot = -(mu*z/(r^3))*(1 + (((Ri/r)^2)*J2*((15/2)*(z/r)^2-(9/2))));

dX = [xDot; yDot; zDot; xDDot; yDDot; zDDot];


end



