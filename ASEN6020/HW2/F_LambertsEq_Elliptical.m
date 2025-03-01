function F = F_LambertsEq_Elliptical(a, mu, s, c, TOF, shortTOF, lt180)
% Definition of Lambert's Equation for elliptical transfers. fsolve 
% iteratively solves this equation for a.
%   Inputs:
%       - a: Current guess of transfer semi-major axis
%       - s: Semi-perimeter of the space triangle for the desired transfer
%       - c: Chord length of the space triangle for the desired transfer
%       - TOF: Desired time of flight for the transfer
%       - shortTOF: Whether the TOF is shorter (1) or longer (0) than
%                   TOFmin
%       - lt180: Whether the desired transfer angle is less than (1) or 
%                greater than (0) 180 degrees
%   Outputs:
%       - F: Function vector for fsolve to iterate
%
%   By: Ian Faber, 10/19/2024
%

n = sqrt(mu/(a^3));
alpha0 = 2*asin(sqrt(s/(2*a)));
beta0 = 2*asin(sqrt((s-c)/(2*a)));

if shortTOF
    alpha = alpha0;
else
    alpha = 2*pi - alpha0;
end

if lt180
    beta = beta0;
else
    beta = -beta0;
end

F = (1/n)*((alpha - beta) - (sin(alpha) - sin(beta))) - TOF;

end