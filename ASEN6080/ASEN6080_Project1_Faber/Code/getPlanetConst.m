function constants = getPlanetConst()
% Function that returns relevant planetary constants for OD problems. 
% Update as necessary!
%   Inputs:
%       - None
%   Outputs:
%       - constants: Structure of constants for use in OD problems.
%                    Made of the following fields:
%                   - mu: G*M for the planet of interest [m^3/s^2]
%                   - J2: Gravitational parameter for body oblateness
%                   - J3: Gravitational parameter for body pear-shapedness
%                   - Ri: Planetary radius, assuming a sphere [m]
%                   - wPlanet: Planet angular velocity vector in ECI 
%                              coordinates [rad/s]
%                   - rho0: Reference atmospheric density [kg/m^3]
%                   - r0: Reference atmospheric radius [m]
%                   - H: Atmospheric altitude cutoff [m]
%
%   By: Ian Faber, 01/30/2025
%

    % Gravitational parameters
constants.mu = 3.986004415e14; % Earth
constants.J2 = 1.082626925638815e-3; % Earth

    % Planetary parameters
constants.Ri = 6378136.3; % Earth
constants.wPlanet = [0; 0; 7.2921158553e-5]; % Earth - 23 hr, 56 min and 4 sec day

    % Atmospheric parameters
constants.rho0 = 3.614e-13;
constants.r0 = 700e3 + constants.Ri;
constants.H = 88667.0;

end
