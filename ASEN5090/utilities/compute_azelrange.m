function [AZ_deg, EL_deg, RANGE] = compute_azelrange(userECEF, satECEF)
% Function that calculates azimuth angle, elevation angle, and range from a
% user's ECEF coordinates to an array of satellite ECEF coordinates
%   Inputs:
%       - userECEF: User ECEF coordinates
%       - satECEF: Array of satellite ECEF coordinates as follows:
%                  satECEF = [sat1_ECEF, sat2_ECEF, ..., satN_ECEF]
%   Outputs:
%       - AZ_deg: Array of azimuth angles to each satellite in satECEF in 
%                 degrees like so:
%                 [AZ_1; AZ_2; ...; AZ_N]
%       - EL_deg: Array of elevation angles to each satellite in satECEF in 
%                 degrees like so:
%                 [EL_1; EL_2; ...; EL_N]
%       - RANGE: Array of ranges to each satellite in satECEF like so:
%                [RANGE_1; RANGE_2; ...; RANGE_N]
%
%   By: Ian Faber, 09/12/2025
%

AZ = [];
EL = [];
RANGE = [];

userLLA = ecef2lla(reshape(userECEF, 1, 3));
C_ECEF2ENU = ECEF2ENU(userLLA(1), userLLA(2));
userENU = C_ECEF2ENU*userECEF;

for k = 1:size(satECEF, 2)
    satECEF_k = satECEF(:,k);
    satENU = C_ECEF2ENU*satECEF_k;

    AZ = [AZ, atan2(satENU(1), satENU(2))];
    EL = [EL, asin(satENU(3)/norm(satENU))];
    RANGE = [RANGE, norm(satENU - userENU)];

end

AZ_deg = rad2deg(AZ);
EL_deg = rad2deg(EL);

end