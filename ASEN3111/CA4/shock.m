function [p2p1, p02p01, M2] = shock(beta12, theta12, M1, gamma)
% shock - Function that performs all necessary pressure and mach 
% calculations for a diamond airfoil oblique shock
%
%   Inputs:
%       beta12: Calculated weak beta angle between state 1 and 2 from the 
%               BTM equation (deg)
%       theta12: Wedge angle between state 1 and 2, i.e. upstream and 
%                downstream (deg)
%       M1: Mach number of state 1, i.e. upstream mach number
%       gamma: Specific heat ratio of the medium the airfoil is flying 
%              through
%
%   Outputs:
%       p2p1: Ratio of downstream to upstream static pressure
%       p02p01: Ratio of downstream to upstream total pressure
%       M2: Downstream mach number
%
%   Author: Ian Faber, 04/19/2023
%   Collaborators: None

% Find normal mach
Mn1 = M1*sind(beta12);

% Calculate output pressures
[~, ~, p2p1, ~, Mn2, p02p01, ~] = flownormalshock(gamma, Mn1, 'mach');

% Calculate output mach
M2 = Mn2/sind(beta12 - theta12);



end