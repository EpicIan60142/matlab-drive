function cubesat = generateCubesat(uMax, axialThrust, startRing, name, marker, color)
% Function that generates CubeSats based on specified parameters
%   Inputs:
%       - uMax: Maximum thrust the CubeSat can produce
%       - axialThrust: Boolean indicating whether to split the specified
%                      max thrust up between 3 axial thrusters: 1 for x,
%                      one for y, and one for z axes
%       - startRing: Starting ring for the CubeSats to sample their initial
%                    state from
%       - name: What to name the CubeSat! Choose something fun ;)
%       - marker: What marker to plot the CubeSat as. Must be a valid
%                 'Marker' LineStyle option
%       - color: What color to plot the CubeSat trajectory as.
%   Outputs:
%       - cubesat: CubeSat output structure with the following fields:
%           - uMax: Max thrust, will either be a 3x1 vector or 1x1 scalar
%                   depending on the value of axialThrust:
%               - axialThrust == true: uMax = [uMax_x; uMax_y; uMax_z]
%               - axialThrust == false: uMax = uMax
%           - X0: Initial state for the current trajectory. At ring 1, 
%                 this is sampled from P0 like so:
%                 X0 = [X_sampled; Y_sampled; Z_sampled; 0; 0; 0]
%           - t0: Initial time for the current trajectory. At ring 1, 
%                 t0 = 0
%           
%
%   By: Ian Faber, 04/19/2025
%
    % Define CubeSat Thrust
if axialThrust
    uMax_i = uMax/3;
    uMax = [uMax_i; uMax_i; uMax_i];
end

    % Define CubeSat starting position
        % Project starting ring into 3-D space
T = startRing.NR; % Get DCM for intertial coordinates from ring coordinates
T = T(:,1:2); % Get projection from nPrime and nDPrime vectors to XYZ space
P = T*startRing.S*T'; % 3sigma Covariance in 3D space

        % Sample starting position
sqrtP = zeros(size(P));
for k = 1:size(P,1)
    sqrtP(k,k) = sqrt(P(k,k));
end
X0_pos = startRing.center + sqrtP*randn(3,1)/3;
X0_vel = zeros(3,1);
X0 = [X0_pos; X0_vel];

    % Assign outputs
cubesat = struct("uMax", uMax, "X0", X0, "t0", 0, "X", [], "t", [], "name", name, "marker", marker, "color", color);

end