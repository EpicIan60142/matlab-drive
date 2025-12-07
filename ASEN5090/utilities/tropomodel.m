function [tropoCorr, idxValid] = tropomodel(zd, elevation, elevationCutoff)
% Function that calculates a simple tropospheric delay model based on
% satellite elevation and an assumed ground station zenith delay
%   Inputs:
%       - zd: Assumed zenith delay for the ground station in meters
%       - elevation: Satellite elevation in radians
%       - elevationCutoff: Elevation angle validity cutoff in degrees, any
%                          elevation angle below this will result in a NaN
%                          result for that elevation
%   Outputs:
%       - tropoCorr: Simple tropospheric correction in meters
%       - idxValid: Vector of valid tropospheric corrections based on the
%                   provided elevation angle cutoff
%
%   By: Ian Faber, 10/03/2025
%

    % Calculate correction
tropoCorr = [];
for k = 1:length(elevation)
    if rad2deg(elevation(k)) > elevationCutoff
        tropoCorr = [tropoCorr; zd/sin(elevation(k))];
    else
        tropoCorr = [tropoCorr; NaN];
    end
end

idxValid = ~isnan(tropoCorr);

end