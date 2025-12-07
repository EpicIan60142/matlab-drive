function C_ECEF2ENU = ECEF2ENU(lat_deg, lon_deg)
% Function that calculates the DCM from ECEF to ENU given latitude and
% longitude in degrees
%   Inputs:
%       - lat_deg: Latitude in degrees North
%       - lon_deg: Longitude in degrees East
%   Outputs:
%       - C_ECEF2ENU: DCM from ECEF to ENU
%
%   By: Ian Faber, 09/12/2025
%

lat_rad = deg2rad(lat_deg);
lon_rad = deg2rad(lon_deg);

C_ECEF2ENU = EA2DCM([pi/2 + lon_rad, pi/2 - lat_rad, 0], [3,1,3]);

end