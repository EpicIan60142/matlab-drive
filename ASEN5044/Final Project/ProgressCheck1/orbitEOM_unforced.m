function dX = orbitEOM_unforced(t,X, mu)
% Nonlinear orbit EOM for ode45 without force inputs and process noise

x = X(1);
xDot = X(2);
y = X(3);
yDot = X(4);

r = sqrt(x^2 + y^2);

xDoubleDot = -mu*x/(r^3);
yDoubleDot = -mu*y/(r^3);

dX = [xDot; xDoubleDot; yDot; yDoubleDot];

end