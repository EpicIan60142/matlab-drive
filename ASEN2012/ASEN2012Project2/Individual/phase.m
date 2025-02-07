function [value, isterminal, direction] = phase(t,X)
% Define events of interest to ODE45, specifically integration termination 
% events
%   Inputs: Time vector, t, and state vector, X, formatted as 
%   [x;z;vx;vz;m;mAir;Vair]
%
%   Outputs: Value to watch, value, whether the value will terminate
%   integration, isterminal, and what direction to watch the value change,
%   direction
%

%Extract current height
z = X(2);

value = z; % which variable to use (hint, this indicates when the z-coordinate hits the GROUND level, not sea level)
isterminal = 1; % terminate integration? (0 or 1)
direction = -1; % test when the value first goes negative (-1) or positive (+1)
end
