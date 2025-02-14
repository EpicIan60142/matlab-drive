function stations = makeStations()
% Function that returns relevant station constants for OD problems. 
% Update as necessary!
%   Outputs:
%       - stations: Structure array of stations with the following
%                   information:
%                   - X0: The initial station position in ECI coordinates
%                   - elMask: The elevation angle mask of the station
%                   - fT: The transmission frequency of the station
%                   - rho: Station range measurements of the spacecraft
%                   - rhoDot: Station range-rate measurements of the
%                             spacecraft
%                   - elAngle: Station elevation angle measurements of the
%                              spacecraft
%                   - tMeas: Vector of times when measurements are made.
%                            There will be one time per measurement in rho,
%                            rhoDot, elAngle, etc.
%                   - sigRho: Range measurement uncertainty in m
%                   - sigRhoDot: Range rate measurement uncertainty in m/s
%                   - R: Measurement uncertainty matrix
%                   - id: Station id as specified in the project statement
%                   - idx: Station index in the overall stations struct
%                   - color: The color to plot this station's measurements
%                            in
%
%   By: Ian Faber, 01/30/2025
%

    % Create empty struct
stations = struct('X0', [], 'elMask', [], 'fT', [], 'rho', [], ...
                  'rhoDot', [], 'elAngle', [], 'tMeas', [], 'sigRho', [], ...
                  'sigRhoDot', [], 'R', {}, 'id', [], 'idx', [], ...
                  'color', []); % Create struct for storing station info

    % Set individual station parameters
elMask = deg2rad([10; 10; 10]);
fT = [8.99e9; 8.99e9; 8.99e9];
sigRho = [0.01; 0.01; 0.01]; % cm -> m
sigRhoDot = [0.001; 0.001; 0.001]; % mm -> m
Xs = [-5127510.0; 3860910.0; 549505.0];
Ys = [-3794160.0; 3238490.0; -1380872.0];
Zs = [0.0; 3898094.0; 6182197.0];
ids = [101; 337; 394];
colors = ['b', 'r', 'k'];

    % Initial spin angle of planet relative to the positive inertial z-axis
theta0 = deg2rad(0);

    % Make initial station information
for k = 1:length(Xs)
    body = [Xs(k); Ys(k); Zs(k)];

    stations(k).X0 = rotZ(theta0)*body;
    stations(k).elMask = elMask(k);
    stations(k).fT = fT(k);
    stations(k).sigRho = sigRho(k);
    stations(k).sigRhoDot = sigRhoDot(k);
    stations(k).id = ids(k);
    stations(k).idx = k;
    stations(k).color = colors(k);
end

end
