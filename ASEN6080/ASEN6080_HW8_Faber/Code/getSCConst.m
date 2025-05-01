function sc = getSCConst()
% Function that returns relevant spacecraft constants for OD problems. 
% Update as necessary!
%   Inputs:
%       - None
%   Outputs:
%       - sc: Structure of spacecraft paremeters at the start of an OD 
%             problem. Made of the following fields:
%           - A: Cross-sectional area [km^2]
%           - m: Mass [kg]
%           - Cd: Drag coefficient [n.d.]
%
%   By: Ian Faber, 02/10/2025
%

    % Physical Parameters
sc.A = 3/(1000^2); % m^2 -> km^2
sc.m = 970;
sc.Cd = 2;

end