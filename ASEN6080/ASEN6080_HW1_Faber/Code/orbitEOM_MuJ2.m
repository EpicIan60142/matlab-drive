function dX = orbitEOM_MuJ2(t,X,mu,Ri)
% Function for the equations of motion for a body in orbit around a central
% body including mu and J2 in the 2 body problem
%   Inputs:
%       - t: Current integration time
%       - X: State vector arranged as follows:
%            [X; Y; Z; Xdot; Ydot; Zdot; J2]
%       - mu: Gravitational parameter for the central body of interest
%       - Ri: Reference radius
%   
%   Outputs:
%       - dX: Rate of change vector for the provided state
%             [Xdot; Ydot; Zdot; Xddot; Yddot; Zddot; J2dot]
%
%   By: Ian Faber, 01/23/2025
%

x = X(1);
y = X(2);
z = X(3);
xDot = X(4);
yDot = X(5);
zDot = X(6);
J2 = X(7);

r = sqrt(x^2 + y^2 + z^2);

xDDot = -(mu*x/(r^3))*(1 + (((Ri/r)^2)*J2*((15/2)*(z/r)^2-(3/2))));
yDDot = -(mu*y/(r^3))*(1 + (((Ri/r)^2)*J2*((15/2)*(z/r)^2-(3/2))));
zDDot = -(mu*z/(r^3))*(1 + (((Ri/r)^2)*J2*((15/2)*(z/r)^2-(9/2))));

J2Dot = 0;

dX = [xDot; yDot; zDot; xDDot; yDDot; zDDot; J2Dot];


end



