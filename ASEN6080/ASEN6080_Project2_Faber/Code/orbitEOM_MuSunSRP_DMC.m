function dX = orbitEOM_MuSunSRP_DMC(t,X,B,u,pConst,scConst)
% Function for the equations of motion for a body in orbit around a central
% body including earth and sun point masses and SRP (cannonball model)
%   Inputs:
%       - t: Current integration time in seconds
%       - X: State vector arranged as follows:
%            [X; Y; Z; Xdot; Ydot; Zdot]
%       - pConst: Planetary constant structure from getPlanetConst.m
%       - scConst: Spacecraft constant structure from getSCConst.m
%   
%   Outputs:
%       - dX: Rate of change vector for the provided state
%             [Xdot; Ydot; Zdot; Xddot; Yddot; Zddot; C_Rdot]
%
%   By: Ian Faber, 04/03/2025
%

    % Extract state
x = X(1);
y = X(2);
z = X(3);
xDot = X(4);
yDot = X(5);
zDot = X(6);
C_R = X(7);
wX = X(8);
wY = X(9);
wZ = X(10);

    % Make vectors
rVec = [x; y; z]; % Position vector from Earth to sc
[RVec_ic, ~, ~, ~] = Ephem(pConst.initEpoch + (t/(24*60*60)), 3, 'EME2000'); % Convert seconds to JD - 1 JD = 24*60*60 sec
RVec_ic = -RVec_ic; % Ephem gives vector from sun to Earth, but want Earth to sun
rVec_isc = RVec_ic - rVec; % Vector from sc to Sun

    % Make constants
r = norm(rVec);
R_ic = norm(RVec_ic);
r_isc = norm(rVec_isc);
Pphi = scConst.Phi/pConst.c;

rSRP = r_isc/149597870.7; % Non-dimensionalize R for SRP model by 1 AU in km

    % Make contributing accelerations
a_pointMassEarth = -(pConst.mu_Earth/(r^3))*rVec;
a_pointMassSun = pConst.mu_Sun*((rVec_isc/(r_isc^3)) - (RVec_ic/(R_ic^3)));
a_SRP = -C_R*(Pphi/(rSRP^2))*scConst.AtoM*(rVec_isc/r_isc);

rDD = a_pointMassEarth + a_pointMassSun + a_SRP + [wX; wY; wZ];

C_Rdot = 0; % C_R is constant

wDot = -B*[wX; wY; wZ];

    % Assign outputs
dX = [xDot; yDot; zDot; rDD; C_Rdot; wDot];


end



