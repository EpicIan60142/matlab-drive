function Y = generateRUDoppler(X, Xs, elMask, fT)
% Function that generates observations of spacecraft from ground stations,
% assuming the only measurements are Range Units and Doppler shift from the 
% DSN. This function also provides the elevation angle, and will only 
% provide a non-NaN measurement when the satellite is visible within the 
% station's elevation mask, which is assumed to be conical for this 
% function.
%   Inputs:
%       - X: Spacecraft state, organized as follows:
%            [x; y; z; xDot; yDot; zDot]
%       - Xs: Ground station state, organized as follows:
%            [x_s; y_s; z_s; xDot_s; yDot_s; zDot_s]
%       - elMask: Ground station elevation mask in RADIANS
%       - fT: DSN Transmit frequency in Hz
%   Outputs:
%       - Y: Measurement vector, organized as follows:
%            [RU; fShift; elAngle]
%
%   By: Ian Faber, 01/26/2025
%

c = 299792; % speed of light in km/s

R = X(1:3);
V = X(4:6);
Rs = Xs(1:3);
Vs = Xs(4:6);

rhoVec = R - Rs;

rho = norm(rhoVec);
rhoDot = dot(R-Rs, V-Vs)/rho;

RU = (221/749)*(rho/c)*fT;
fShift = -2*(rhoDot/c)*fT;

elAngle = asin(dot(rhoVec, Rs)/(rho*norm(Rs)));

if abs(elAngle) >= elMask && elAngle > 0
    Y = [RU; fShift; elAngle];
else
    Y = [NaN; NaN; NaN];
end

end