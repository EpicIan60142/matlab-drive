function F = F_LambertsEq_Hyperbolic(a, mu, s, c, TOF, shortTOF, lt180)
% Definition of Lambert's Equation for hyperbolic transfers. fsolve 
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

n = sqrt(mu/(abs(a)^3));
alpha0 = 2*asinh(sqrt(s/(2*abs(a))));
beta0 = 2*asinh(sqrt((s-c)/(2*abs(a))));

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

% Only one of these function definitions is ever valid during one run of
% Lambert's problem!
if lt180 % Transfer angle < 180 degrees
    F = (1/n)*(sinh(alpha) - alpha - (sinh(beta) - beta)) - TOF;
else
    F = (1/n)*(sinh(alpha) - alpha + (sinh(beta) - beta)) - TOF;
end

end