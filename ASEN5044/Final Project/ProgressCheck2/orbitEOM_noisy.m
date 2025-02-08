function dX = orbitEOM_noisy(t,X, mu, w)
% Nonlinear orbit EOM for ode45 with process noise

x = X(1);
xDot = X(2);
y = X(3);
yDot = X(4);

r = sqrt(x^2 + y^2);

xDoubleDot = -mu*x/(r^3) + w(1);
yDoubleDot = -mu*y/(r^3) + w(2);

dX = [xDot; xDoubleDot; yDot; yDoubleDot];

end