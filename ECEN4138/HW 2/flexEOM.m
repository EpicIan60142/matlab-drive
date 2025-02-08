function dX = flexEOM(t, X, const)
% EOM function for simulating a simplified 1-D flexible system with ode45
%   Inputs:
%       t: time [sec]
%       X: state vector
%           [ x; y; vx; vy ]
%       const: vector of constants for simulation
%           [M; m; b1; b2; k; u] -> If not simulating air drag, b2 = 0
%
%   Outputs:
%       dX: rate of change vector
%           [ vx; vy; ax; ay ]
%
%   By: Ian Faber, 09/11/2023
%

M = const(1);
m = const(2);
b1 = const(3);
b2 = const(4);
k = const(5);
u = const(6);

x = X(1);
y = X(2);
vx = X(3);
vy = X(4);

% if t > 30
%     u = 0;
% end

ax = (1/m)*(b1*vy + k*y - b1*vx - k*x);
ay = (1/M)*(u + b1*vx + k*x - (b1+b2)*vy - k*y );

dX = [vx; vy; ax; ay];

end