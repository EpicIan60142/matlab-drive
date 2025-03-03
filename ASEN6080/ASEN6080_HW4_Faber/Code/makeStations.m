function stations = makeStations(pConst)
% Function that returns relevant station constants for OD problems. 
% Update as necessary!
%   Inputs:
%       - pConst: Planetary constant structure, generated with
%                 getPlanetConst.m
%   Outputs:
%       - stations: Structure array of stations with the following
%                   information:
%                   - X0: The initial station position in ECI coordinates
%                   - elMask: The elevation angle mask of the station
%                   - fT: The transmission frequency of the station
%                   - Xs: The state of the station
%                   - rho: Station range measurements of the spacecraft
%                   - rhoDot: Station range-rate measurements of the
%                             spacecraft
%                   - elAngle: Station elevation angle measurements of the
%                              spacecraft
%                   - RU: Station range unit measurements of the spacecraft
%                   - fShift: Station doppler shift measurements of the
%                             spacecraft
%                   - tMeas: Vector of times when measurements are made.
%                            There will be one time per measurement in rho,
%                            rhoDot, elAngle, etc.
%                   - sigRho: Range measurement uncertainty in km
%                   - sigRhoDot: Range rate measurement uncertainty in km/s
%                   - R: Measurement uncertainty matrix
%                   - id: Station id in order of creation
%                   - color: The color to plot this station's measurements
%                            in
%
%   By: Ian Faber, 01/30/2025
%

    % Create empty struct
stations = struct('X0', [], 'elMask', [], 'fT', [], 'Xs', [], 'rho', [], ...
                  'rhoDot', [], 'elAngle', [], 'RU', [], 'fShift', [], ...
                  'tMeas', [], 'sigRho', [], 'sigRhoDot', [], 'R', {}, ...
                  'id', [], 'color', []); % Create struct for storing station info

    % Set individual station parameters
elMask = deg2rad([10; 10; 10]);
fT = [8.99e9; 8.99e9; 8.99e9];
sigRho = [1e-3; 1e-3; 1e-3];
sigRhoDot = [1e-6; 1e-6; 1e-6];
latitudes = deg2rad([-35.39833; 40.42722; 35.247164]); % phi
longitudes = deg2rad([148.981944; 355.749444; 243.205]); % lambda
colors = ['b', 'r', 'k'];

    % Initial spin angle of planet relative to the positive inertial z-axis
theta0 = deg2rad(122);

    % Make initial station information
for k = 1:length(latitudes)
    spherical = [pConst.Ri; latitudes(k); longitudes(k)];
    body = toBodyFromSpherical(spherical);

    stations(k).X0 = rotZ(theta0)*body;
    stations(k).elMask = elMask(k);
    stations(k).fT = fT(k);
    stations(k).sigRho = sigRho(k);
    stations(k).sigRhoDot = sigRhoDot(k);
    stations(k).id = k;
    stations(k).color = colors(k);
end

end
