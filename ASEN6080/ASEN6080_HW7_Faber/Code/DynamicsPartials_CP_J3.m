function B = DynamicsPartials_CP_J3(X, pConst)
% Function that outputs the acceleration partials matrix for orbital 
% dynamics including contributions from Mu, J2 with respect to the consider
% parameter, J3
%   Inputs:
%       - X: System state in SI units 
%           [(x,y,z) -> km, (xDot, yDot, zDot) -> km/s
%           X = [x; y; z; xDot; yDot; zDot]
%       - pConst: Planetary constants structure as defined in
%                 getPlanetConst.m
%
%   Outputs:
%       - B: Acceleration partials matrix with respect to consider
%            parameter J3
%
%   By: Ian Faber, 04/06/2025
%

    % Extract states
x = X(1);
y = X(2);
z = X(3);

    % Extract constants
Ri = pConst.Ri;
mu = pConst.mu;

    % Define range
r = sqrt(x^2 + y^2 + z^2);

    % Define non-trivial partials
delXddDelJ3 = (-mu*x/r^3)*((Ri/r)^3)*((35/2)*(z/r)^3 - (15/2)*(z/r));

delYddDelJ3 = (-mu*y/r^3)*((Ri/r)^3)*((35/2)*(z/r)^3 - (15/2)*(z/r));

delZddDelJ3 = -(mu/r^2)*((Ri/r)^3)*((35/2)*(z/r)^4 - 15*(z/r)^2 + (3/2));

    % Create matrix blocks for convenience
constantsBlock = [
                    delXddDelJ3;
                    delYddDelJ3;
                    delZddDelJ3
                 ];

    % Construct B matrix
B = [
        zeros(3,1); % delFdelJ3 is 0 wrt xDot, yDot, zDot
        constantsBlock
    ];

end




