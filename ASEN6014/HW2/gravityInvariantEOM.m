function dX = gravityInvariantEOM(t, X, const)
% Function to be passed into ode45 to calculate the evolution of the
% invariants of motion for a spring mass damper with gravity disturbance
% acceleration
%   Inputs:
%       - t: Current simulation time in sec
%       - X: State vector composed of [e1; e2]
%       - const: Constants structure with fields m, g, and k
%   Outputs:
%       - dX: Rate of change vector [de1; de2]
%
%   By: Ian Faber, 09/13/2025
%

% Extract state
e1 = X(1);
e2 = X(2);

% Extract constants
m = const.m;
g = const.g;
k = const.k;

% Calculate angular rate
w = sqrt(k/m);

% Calculate rate of change
dX = [(g/w)*sin(w*t); -g*cos(w*t)];


end