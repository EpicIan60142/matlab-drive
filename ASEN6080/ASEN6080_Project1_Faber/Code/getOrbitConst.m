function orbital = getOrbitConst()
% Function that returns relevant planetary constants for OD problems. 
% Update as necessary!
%   Inputs:
%       - None
%   Outputs:
%       - constants: Structure of orbital elements at the start of an OD 
%                    problem. Made of the following fields:
%                   - a: Semi-major axis [km]
%                   - e: Eccentricity
%                   - i: Inclination [rad]
%                   - RAAN: Right ascension of the ascending node [rad]
%                   - argPeri: Argument of periapsis [rad]
%                   - truAnom: True anomaly [rad]
%
%   By: Ian Faber, 01/30/2025
%

orbital.a = 800+6378;
orbital.e = deg2rad(0); % unknown
orbital.i = deg2rad(98.6);
orbital.RAAN = deg2rad(0); % unknown
orbital.argPeri = deg2rad(0); % unknown
orbital.truAnom = deg2rad(0); % unknown

end
