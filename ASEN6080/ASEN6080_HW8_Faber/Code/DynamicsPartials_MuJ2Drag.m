function A = DynamicsPartials_MuJ2Drag(X, pConst, scConst)
% Function that outputs the acceleration partials matrix for orbital 
% dynamics including contributions from mu, J2, and drag.
%   Inputs:
%       - X: System state in km
%           X = [x; y; z; xDot; yDot; zDot]
%       - pConst: Planetary constants structure as defined by
%                 getPlanetConst.m
%
%   Outputs:
%       - A: Acceleration partials matrix
%
%   By: Ian Faber, 01/22/2025
%

    % Extract states
x = X(1);
y = X(2);
z = X(3);
xDot = X(4);
yDot = X(5);
zDot = X(6);

    % Extract constants
mu = pConst.mu;
J2 = pConst.J2;
Ri = pConst.Ri;
rho0 = pConst.rho0;
r0 = pConst.r0;
H = pConst.H;

Asc = scConst.A;
m = scConst.m;
Cd = scConst.Cd;

    % Define range and speed relative to ECI frame
r = sqrt(x^2 + y^2 + z^2);
v = sqrt(xDot^2 + yDot^2 + zDot^2);

    % Calculate atmospheric density
rho = calcDensity(pConst, r);

    % Define helper constants
dragConst = 0.5*(Cd*Asc/m);
atmosConst = (rho0/(H*r))*exp(-(r-r0)/H);

    % Define non-trivial partials
delXddDelX = -mu*( ((r^2 - 3*x^2)/r^5) + ((Ri^2*J2)/2)*(((15*(x^2 + z^2))/r^7) - ((105*x^2*z^2)/r^9) - (3/r^5)) ) + dragConst*v*xDot*x*atmosConst;
delXddDelY = mu*( ((3*x*y)/r^5) - ((Ri^2*J2*x*y)/2)*((15/r^7) - ((105*z^2)/r^9)) ) + dragConst*v*xDot*y*atmosConst;
delXddDelZ = mu*( ((3*x*z)/r^5) - ((Ri^2*J2*x*z)/2)*((45/r^7) - ((105*z^2)/r^9)) ) + dragConst*v*xDot*z*atmosConst;

delXddDelXd = -dragConst*rho*(((xDot^2)/v) + v);
delXddDelYd = -dragConst*rho*xDot*yDot/v;
delXddDelZd = -dragConst*rho*xDot*zDot/v;

delYddDelX = mu*( ((3*x*y)/r^5) - ((Ri^2*J2*x*y)/2)*((15/r^7) - ((105*z^2)/r^9)) ) + dragConst*v*yDot*x*atmosConst;
delYddDelY = -mu*( ((r^2 - 3*y^2)/r^5) + ((Ri^2*J2)/2)*(((15*(y^2 + z^2))/r^7) - ((105*y^2*z^2)/r^9) - (3/r^5)) ) + dragConst*v*yDot*y*atmosConst;
delYddDelZ = mu*( ((3*y*z)/r^5) - ((Ri^2*J2*y*z)/2)*((45/r^7) - ((105*z^2)/r^9)) ) + dragConst*v*yDot*z*atmosConst;

delYddDelXd = -dragConst*rho*yDot*xDot/v;
delYddDelYd = -dragConst*rho*(((yDot^2)/v) + v);
delYddDelZd = -dragConst*rho*yDot*zDot/v;

delZddDelX = mu*( ((3*x*z)/r^5) - ((Ri^2*J2*x*z)/2)*((45/r^7) - ((105*z^2)/r^9)) ) + dragConst*v*zDot*x*atmosConst;
delZddDelY = mu*( ((3*y*z)/r^5) - ((Ri^2*J2*y*z)/2)*((45/r^7) - ((105*z^2)/r^9)) ) + dragConst*v*zDot*y*atmosConst;
delZddDelZ = -mu*( ((r^2 - 3*z^2)/r^5) + ((Ri^2*J2)/2)*(((90*z^2)/r^7) - ((105*z^4)/r^9) - (9/r^5)) ) + dragConst*v*zDot*z*atmosConst;

delZddDelXd = -dragConst*rho*zDot*xDot/v;
delZddDelYd = -dragConst*rho*zDot*yDot/v;
delZddDelZd = -dragConst*rho*(((zDot^2)/v) + v);

delXddDelMu = (-x/r^3)*(1 + ((Ri/r)^2)*J2*((15/2)*(z/r)^2 - (3/2)));
delYddDelMu = (-y/r^3)*(1 + ((Ri/r)^2)*J2*((15/2)*(z/r)^2 - (3/2)));
delZddDelMu = (-z/r^3)*(1 + ((Ri/r)^2)*J2*((15/2)*(z/r)^2 - (9/2)));

delXddDelJ2 = (-mu*x/r^3)*((Ri/r)^2)*((15/2)*(z/r)^2 - (3/2));
delYddDelJ2 = (-mu*y/r^3)*((Ri/r)^2)*((15/2)*(z/r)^2 - (3/2));
delZddDelJ2 = (-mu*z/r^3)*((Ri/r)^2)*((15/2)*(z/r)^2 - (9/2));

delXddDelCd = -0.5*rho*(Asc/m)*v*xDot;
delYddDelCd = -0.5*rho*(Asc/m)*v*yDot;
delZddDelCd = -0.5*rho*(Asc/m)*v*zDot;

    % Create matrix blocks for convenience
dynamicsBlock_1 = [
                    delXddDelX, delXddDelY, delXddDelZ;
                    delYddDelX, delYddDelY, delYddDelZ;
                    delZddDelX, delZddDelY, delZddDelZ
                  ];

dynamicsBlock_2 = [
                    delXddDelXd, delXddDelYd, delXddDelZd;
                    delYddDelXd, delYddDelYd, delYddDelZd;
                    delZddDelXd, delZddDelYd, delZddDelZd
                  ];

% constantsBlock = [
%                     delXddDelMu, delXddDelJ2, delXddDelCd;
%                     delYddDelMu, delYddDelJ2, delYddDelCd;
%                     delZddDelMu, delZddDelJ2, delZddDelCd
%                  ];

% stationsBlock = [];
% for k = 1:3
%     stationsBlock = blkdiag(stationsBlock, tilde(pConst.wPlanet));
% end

    % Construct A matrix
A = [
        zeros(3,3), eye(3);
        dynamicsBlock_1, dynamicsBlock_2;
    ];

end




