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
%                   - Xs: The state of the station
%                   - rho: Station range measurements of the spacecraft
%                   - rhoDot: Station range-rate measurements of the
%                             spacecraft
%                   - elAngle: Station elevation angle measurements of the
%                              spacecraft
%                   - tMeas: Vector of times when measurements are made.
%                            There will be one time per measurement in rho,
%                            rhoDot, elAngle, etc.
%                   - sigRho: Range measurement uncertainty in km
%                   - sigRhoDot: Range rate measurement uncertainty in km/s
%                   - R: Measurement uncertainty matrix
%                   - id: Station id for labeling, string
%                   - idx: Station index in order of creation
%                   - color: The color to plot this station's measurements
%                            in
%
%   By: Ian Faber, 01/30/2025
%

    % Create empty struct
stations = struct('X0', [], 'elMask', [], 'Xs', [], 'rho', [], 'rhoDot', [], ...
                  'elAngle', [],  'tMeas', [], 'sigRho', [], 'sigRhoDot', [], ...
                  'R', {}, 'id', [], 'idx', [], 'color', []); % Create struct for storing station info

    % Set individual station parameters
elMask = deg2rad([10; 10; 10]);
sigRho = [5e-3; 5e-3; 5e-3]; % 5 m uncertainty in range, km
sigRhoDot = [0.5e-6; 0.5e-6; 0.5e-6]; % 0.5 mm/s uncertainty in range rate, km/s
latitudes = deg2rad([-35.398333; 40.427222; 35.247164]); % phi
longitudes = deg2rad([148.981944; 355.749444; 243.205]); % lambda
altitudes = [0.691750; 0.834539; 1.07114904]; % km
ids = ['DSS 34', 'DSS 65', 'DSS 13'];
colors = ['b', 'r', 'k'];

    % Initial spin angle of planet relative to the positive inertial z-axis
theta0 = deg2rad(0);

    % Make initial station information
for k = 1:length(latitudes)
    spherical = [pConst.Ri + altitudes(k); latitudes(k); longitudes(k)];
    body = toBodyFromSpherical(spherical);

    stations(k).X0 = rotZ(theta0)*body;
    stations(k).elMask = elMask(k);
    stations(k).sigRho = sigRho(k);
    stations(k).sigRhoDot = sigRhoDot(k);
    stations(k).id = ids(k);
    stations(k).idx = k;
    stations(k).color = colors(k);
end

end
