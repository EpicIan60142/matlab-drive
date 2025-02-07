function [xAirfoil, yAirfoil, xVort, yVort, dydxTheta] = NACA_Airfoils(m, p, t, c, N)
%   NACA_Airfoils - 4-digit NACA Airfoil Panel Calculator Function
%       Calculates panels for any 4-digit NACA airfoil
%       Inputs:
%           m: Maximum camber (% x/c)
%           p: Location of max camber (10 * % x/c)
%           t: Thickness (% x/c)
%           c: Chord length (m)
%           N: Number of panels
%       Outputs:
%           xAirfoil: Vector containing the x location of airfoil boundary 
%                     points structured: [xUpper; xLower]
%           yAirfoil: Vector containing the y location of airfoil boundary 
%                     points structured: [yUpper; yLower]
%
%       Example: To calculate 100 panels for a NACA 4315 airfoil with a
%                chord length of 10 m, call the function like so:
%
%           [x, y] = NACA_Airfoils(4,3,15,10,100)
%
%       Author: Ian Faber, 03/12/2023

% Convert NACA numbers
m = m/100;
p = p/10;
t = t/100;

% Define xi mapping
xi = @(dydx) atan(dydx);

% Define base x vector
x = linspace(0,c,floor((N/2)+1));

a = ((x >= 0) & (x < p*c)); % Before max camber
b = ((x >= p*c) & (x <= c)); % After max camber

% Define thickness distribution
yt = (t/0.2)*c*(0.2969*sqrt(x/c) - 0.126*(x/c) - 0.3516*(x/c).^2 + 0.2843*(x/c).^3 -0.1036*(x/c).^4);

% Define mean camber line
if p ~= 0
    yc = a.*(m*(x/p^2).*(2*p - x/c)) + b.*(m*((c-x)/(1-p)^2).*(1 + x/c - 2*p));
else
    yc = (m*((c-x)/(1-p)^2).*(1 + x/c - 2*p));
end

% Define derivative of the camber line
if p ~= 0
    dydx = a*(m/p^2).*(2*p - (2*x)/c) + b*(m/(1-p)^2).*(2*p - (2*x)/c);
else
    dydx = (m/(1-p)^2).*(2*p - (2*x)/c);
end

% Define mapped derivative of camber line
if p ~= 0
    dydxTheta = @(theta) a.*(m/p^2).*(2*p - 2*(0.5*(1-cos(theta)))) + b.*(m/(1-p)^2).*(2*p - 2*(0.5*(1-cos(theta))));
else
    dydxTheta = @(theta) (m/(1-p)^2).*(2*p - 2*(0.5*(1-cos(theta))));
end

% Calculate coordinates
xUpper = x - yt.*sin(xi(dydx));
xLower = x + yt.*sin(xi(dydx));

yUpper = yc + yt.*cos(xi(dydx));
yLower = yc - yt.*cos(xi(dydx));

% Compile normal coordinates
xAirfoil = [xUpper; xLower];
yAirfoil = [yUpper; yLower];

% Compile wrapped coordinates for vortex panel method
xVort = [xLower(end:-1:1), xUpper(2:end)];
yVort = [yLower(end:-1:1), yUpper(2:end)];

end
