function constants = getPlanetConst()
% Function that returns relevant planetary constants for OD problems. 
% Update as necessary!
%   Inputs:
%       - None
%   Outputs:
%       - constants: Structure of constants for use in OD problems.
%                    Made of the following fields:
%                   - mu_Earth: G*M for Earth [km^3/s^2]
%                   - mu_Sun: G*M for the Sun [km^3/s^2]
%                   - Ri: Planetary radius of Earth, assuming a sphere [km]
%                   - wEarth: Earth angular velocity about the positive 
%                             inertial z-axis [rad/s]
%                   - initEpoch: Initial epoch of the ephemerides being
%                                used for Earth [JD]
%                   - c: Speed of light [km/s]
%
%   By: Ian Faber, 04/03/2025
%
    % Time parameters
constants.initEpoch = 2456296.25; % JD

    % Gravitational parameters
[~,~,mu_Earth,mu_Sun] = Ephem(constants.initEpoch, 3, 'EME2000');

constants.mu_Earth = mu_Earth;
constants.mu_Sun = mu_Sun;

    % Planetary parameters
constants.Ri = 6378.1363; % Earth, km
constants.RSOI = 925000; % Earth, km
constants.wEarth = 7.29211585275553e-5; % Earth, rad/s

    % Other parameters
constants.c = 299792.458; % Speed of light, km/s

end
