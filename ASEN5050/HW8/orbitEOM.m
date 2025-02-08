function dX = orbitEOM(t, X, mu)
% Orbit EOM function for ODE45
%   - Ian Faber, 11/20/2024
%
x = X(1);
y = X(2);
z = X(3);
xDot = X(4);
yDot = X(5);
zDot = X(6);

R = [x; y; z];
r = norm(R);
accel = -mu*R/r^3;

dX = [xDot; yDot; zDot; accel];

end