function [value, isterminal, direction] = SOICheck(t,X,pConst)
% Define events of interest to ODE45, specifically integration termination 
% events
%   Inputs: Time vector, t, and state vector, X, formatted as 
%           [x; y; z; vx; vy; vz; C_R]
%
%   Outputs: Value to watch, value, whether the value will terminate
%   integration, isterminal, and what direction to watch the value change,
%   direction
%

x = X(1);
y = X(2);
z = X(3);

r = sqrt(x^2 + y^2 + z^2);

checkSOI = r - 3*pConst.RSOI; % Stop integration when we hit 3*R_SOI for Earth

value = checkSOI; % Check for this to hit 0
isterminal = 1; % Stop integration if condition is met
direction = -1; % Condition met when r - 3*RSOI goes negative

end