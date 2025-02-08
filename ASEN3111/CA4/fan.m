function [p2p02, M2] = fan(theta12, M1, gamma)
% fan - Function that performs all necessary pressure, mach, and angle 
% calculations for a diamond airfoil Prandtl-Meyer expansion fan
%
%   Inputs:
%       theta12: Expansion angle between states 1 and 2 (deg)
%       M1: Upstream mach number
%       gamma: Specific heat ratio for the medium the airfoil is flying
%       through
%
%   Outputs:
%       p2p02: Ratio of downstream static to total pressure
%       M2: Downstream mach number
%
%   Author: Ian Faber, 04/19/2023
%   Collaborators: None

% Calculate nu's
[~, nu1, ~] = flowprandtlmeyer(gamma, M1, 'mach');
nu2 = nu1 + theta12;

% Calculate outputs
[M2, ~, ~] = flowprandtlmeyer(gamma, nu2, 'nu');
[~, ~, p2p02, ~, ~] = flowisentropic(gamma, M2, 'mach');

end