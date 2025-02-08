function [ts] = calcImpact(t0,vy0,y0) 

% g: gravitational accerelation as local variable
g = 9.81;

% Rearrange for the standard form
a = -0.5*g;
b = g*t0 + vy0;
c = -0.5*g*t0^2 - vy0*t0 + y0;

% Solve the quadratic equation (use the negative root)
ts = (-b - sqrt(b^2 - 4*a*c))/(2*a);

end