function Y = generateRngRngRate(X, statVis, elMask, pConst, ignoreMask)
% Function that generates observations of spacecraft from ground stations,
% assuming the only measurements are simplified range and range rate. This
% function also provides the elevation angle, and will only provide a
% non-NaN measurement when the satellite is visible within the station's
% elevation mask, which is assumed to be conical for this function.
%   Inputs:
%       - X: Spacecraft state, organized as follows:
%            [x; y; z; xDot; yDot; zDot; mu; J2; Cd; statPos_1; statPos_2; 
%             statPos_3]
%       - statVis: Visible station index that is generating the measurement
%       - elMask: Ground station elevation mask in RADIANS
%       - pConst: Planetary constants structure as defined by
%                 getPlanetConst.m
%       - ignoreMask: Boolean indicating if the elevation mask is ignored
%                     (1) or enforced (0). If not provided, defaults to 0.
%   Outputs:
%       - Y: Measurement vector, organized as follows:
%            [rho; rhoDot; elAngle]
%
%   By: Ian Faber, 02/11/2025
%

if ~exist("ignoreMask", 'var')
    ignoreMask = false;
end

R = X(1:3);
V = X(4:6);

idx = statVis;
Rs = X((7+3*idx):(7+3*idx+2));
Vs = cross(pConst.wPlanet, Rs);

rhoVec = R - Rs;

rho = norm(rhoVec);
rhoDot = dot(R-Rs, V-Vs)/rho;

elAngle = asin(dot(rhoVec, Rs)/(rho*norm(Rs)));

if (abs(elAngle) >= elMask && elAngle > 0) || ignoreMask
    Y = [rho; rhoDot; elAngle];
else
    Y = [NaN; NaN; NaN];
end

end