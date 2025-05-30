function constants = getPlanetConst()
% Function that returns relevant planetary constants for OD problems. 
% Update as necessary!
%   Inputs:
%       - None
%   Outputs:
%       - constants: Structure of constants for use in OD problems.
%                    Made of the following fields:
%                   - mu: G*M for the planet of interest [km^3/s^2]
%                   - J2: Gravitational parameter for body oblateness
%                   - Ri: Planetary radius, assuming a sphere [km]
%                   - wPlanet: Planet angular velocity about the positive 
%                              inertial z-axis [rad/s]
%                   - rho0: Reference atmospheric density [kg/km^3]
%                   - r0: Reference atmospheric radius [km]
%                   - H: Atmospheric altitude cutoff [km]
%
%   By: Ian Faber, 01/30/2025
%

    % Gravitational parameters
constants.mu = 398600.4415; % Earth
constants.J2 = 1.08264e-3; % Earth

    % Planetary parameters
constants.Ri = 6378; % Earth
constants.wPlanet = (2*pi)/(24*60*60); % Earth

    % Atmospheric parameters
constants.rho0 = 3.614e-13*(1000^3); % kg/m^3 -> kg/km^3
constants.r0 = 700 + constants.Ri;
constants.H = 88667.0/1000; % m -> km

end
