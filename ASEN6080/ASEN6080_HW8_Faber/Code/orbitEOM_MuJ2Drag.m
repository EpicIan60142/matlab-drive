function dX = orbitEOM_MuJ2Drag(t,X,pConst,scConst)
% Function for the equations of motion for a body in orbit around a central
% body including mu, J2, and drag in the 2 body problem
%   Inputs:
%       - t: Current integration time
%       - X: State vector arranged as follows:
%            [X; Y; Z; Xdot; Ydot; Zdot]
%       - pConst: Planetary constants structure as defined by
%                 getPlanetConst.m
%       - scConst: Spacecraft constants structure as defined by
%                  getSCConst.m
%   
%   Outputs:
%       - dX: Rate of change vector for the provided state
%
%   By: Ian Faber, 02/11/2025
%
    % Extract state
x = X(1);
y = X(2);
z = X(3);
xDot = X(4);
yDot = X(5);
zDot = X(6);

    % Extract constants
Cd = scConst.Cd;
Asc = scConst.A;
m = scConst.m;
dragConst = Cd*Asc/m;
Ri = pConst.Ri;
mu = pConst.mu;
J2 = pConst.J2;

    % Calculate itermediate values
r = sqrt(x^2 + y^2 + z^2);
v = sqrt(xDot^2 + yDot^2 + zDot^2);
rho = calcDensity(pConst, r);

    % Calculate accelerations
xDDot = -(mu*x/(r^3))*(1 + (((Ri/r)^2)*J2*((15/2)*(z/r)^2-(3/2)))) - 0.5*rho*dragConst*v*xDot;
yDDot = -(mu*y/(r^3))*(1 + (((Ri/r)^2)*J2*((15/2)*(z/r)^2-(3/2)))) - 0.5*rho*dragConst*v*yDot;
zDDot = -(mu*z/(r^3))*(1 + (((Ri/r)^2)*J2*((15/2)*(z/r)^2-(9/2)))) - 0.5*rho*dragConst*v*zDot;

    % Assign rates of change
dX = [xDot; yDot; zDot; xDDot; yDDot; zDDot];

end



