function cartesian = convOrbitalToCartesian(orbital)
% Function that converts orbital elements into cartesian states
%   Inputs:
%       - orbital: Vector of orbital elements as follows, where all angles 
%                  are in radians and a and mu are in comparable units 
%                  (i.e. km and km^3/s^2):
%                  [mu; a; e; i; RAAN; argPeri; truAnom]
%   Outputs:
%       - cartesian: Vector of cartesian states as follows, in the same 
%                    units as a and mu:
%                    [x; y; z; xDot; yDot; zDot]
%
%   By: Ian Faber, 01/23/2025
%

mu = orbital(1);
a = orbital(2);
e = orbital(3);
i = orbital(4);
RAAN = orbital(5);
argPeri = orbital(6);
truAnom = orbital(7);

r = (a*(1-e^2))/(1+e*cos(truAnom));
p = a*(1-e^2);
h = sqrt(mu*a*(1-e^2));

vr = (h*e/p)*sin(truAnom);
vtheta = h/r;

Q11 = cos(argPeri)*cos(RAAN) - sin(argPeri)*sin(RAAN)*cos(i);
Q12 = -sin(argPeri)*cos(RAAN) - cos(argPeri)*sin(RAAN)*cos(i);
Q21 = cos(argPeri)*sin(RAAN) + sin(argPeri)*cos(RAAN)*cos(i);
Q22 = -sin(argPeri)*sin(RAAN) + cos(argPeri)*cos(RAAN)*cos(i);
Q31 = sin(argPeri)*sin(i);
Q32 = cos(argPeri)*sin(i);

xDotStar = vr*cos(truAnom) - vtheta*sin(truAnom);
yDotStar = vr*sin(truAnom) + vtheta*cos(truAnom);

x = r*cos(truAnom)*Q11 + r*sin(truAnom)*Q12;
y = r*cos(truAnom)*Q21 + r*sin(truAnom)*Q22;
z = r*cos(truAnom)*Q31 + r*sin(truAnom)*Q32;
xDot = xDotStar*Q11 + yDotStar*Q12;
yDot = xDotStar*Q21 + yDotStar*Q22;
zDot = xDotStar*Q31 + yDotStar*Q32;

cartesian = [x; y; z; xDot; yDot; zDot];

end