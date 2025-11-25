function dX = springInvariantEOM(t, X, const)
% Function to be passed into ode45 to calculate the evolution of the
% invariants of motion for a spring mass damper with cubic disturbance
% accelration
%   Inputs:
%       - t: Current simulation time in sec
%       - X: State vector composed of [e1; e2]
%       - const: Constants structure with fields k1 and k3
%   Outputs:
%       - dX: Rate of change vector [de1; de2]
%
%   By: Ian Faber, 09/13/2025
%

% Extract state
e1 = X(1);
e2 = X(2);

% Extract constants
k1 = const.k1;
k3 = const.k3;

% Calculate unperturbed x at this time
x = e1*cos(sqrt(k1)*t) + (e2/sqrt(k1))*sin(sqrt(k1)*t);

% Disturbance acceleration
ad = -k3*x^3;

% Calculate rate of change
dX = [-ad*((1/sqrt(k1))*sin(sqrt(k1)*t)); ad*cos(sqrt(k1)*t)];


end