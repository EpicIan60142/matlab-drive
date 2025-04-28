function Aprime = DynamicsPartials_MuSunSRP_DMC(X, t, B, pConst, scConst)
% Function that outputs the acceleration partials matrix for orbital 
% dynamics including 3rd body perturbation from the Sun and SRP
%   Inputs:
%       - X: System state in SI units 
%           [(x,y,z) -> km, (xDot, yDot, zDot) -> km/s
%           X = [x; y; z; xDot; yDot; zDot; C_R]
%       - t: Current integration time since the epoch in sec
%       - pConst: Planetary constants structure as defined by
%                 getPlanetConst.m
%       - scConst: Spacecraft constants structure as defined by
%                  getSCConst.m
%
%   Outputs:
%       - A: Acceleration partials matrix
%
%   By: Ian Faber, 04/14/2025
%

    % Extract states
x = X(1);
y = X(2);
z = X(3);
% xDot = X(4);
% yDot = X(5);
% zDot = X(6);
C_R = X(7);

    % Make vectors
rVec = [x; y; z]; % Position vector from Earth to sc
[RVec_ic, ~, ~, ~] = Ephem(pConst.initEpoch + (t/(24*60*60)), 3, 'EME2000'); % Convert seconds to JD - 1 JD = 24*60*60 sec
RVec_ic = -RVec_ic; % Ephem gives vector from sun to Earth, but want Earth to sun
rVec_isc = RVec_ic - rVec; % Vector from sc to Sun

    % Make constants
r = norm(rVec);
% R_ic = norm(RVec_ic);
r_isc = norm(rVec_isc);
Pphi = scConst.Phi/pConst.c;

d = 149597870.7; % 1 AU in km

    % Define unit vectors
rHat = rVec/r;
r_iscHat = rVec_isc/r_isc;

    % Define non-trivial partials
pointMassEarth = (-pConst.mu_Earth/r^3)*(eye(3) - 3*(rHat*rHat'));
pointMassSun = (pConst.mu_Sun/r_isc^3)*(3*(r_iscHat*r_iscHat') - eye(3));
SRP = ((-C_R*d^2*Pphi*scConst.AtoM)/r_isc^3)*(3*(r_iscHat*r_iscHat') - eye(3));

    % Create matrix blocks for convenience
dynamicsBlock = pointMassEarth + pointMassSun + SRP;

constantsBlock = -d^2*Pphi*scConst.AtoM*(rVec_isc/r_isc^3);

    % Construct A matrix
A = [
        zeros(3,3), eye(3), zeros(3,1);
        dynamicsBlock, zeros(3,3), constantsBlock;
        zeros(1,7)
    ];

    % Construct D matrix
D = [zeros(3,3); eye(3); zeros(1,3)];

    % Construct Aprime
Aprime = [
            A, D;
            zeros(3,7), -B
         ];

end




