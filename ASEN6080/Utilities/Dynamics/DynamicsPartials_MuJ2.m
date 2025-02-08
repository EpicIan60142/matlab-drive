function A = DynamicsPartials_MuJ2(X, mu, J2, Ri)
% Function that outputs the acceleration partials matrix for orbital 
% dynamics including contributions from mu and J2.
%   Inputs:
%       - X: System state in SI units 
%           [(x,y,z) -> km, (xDot, yDot, zDot) -> km/s
%           X = [x; y; z; xDot; yDot; zDot]
%       - mu: Gravitational parameter for planetary point mass [km^3/s^2]
%       - J2: Gravitational parameter for planetary oblateness
%       - Ri: Reference radius [km]
%
%   Outputs:
%       - A: Acceleration partials matrix
%
%   By: Ian Faber, 01/22/2025
%

    % Extract states
x = X(1);
y = X(2);
z = X(3);
% xDot = X(4);
% yDot = X(5);
% zDot = X(6);

    % Define range
r = sqrt(x^2 + y^2 + z^2);

    % Define non-trivial partials
delXddDelX = -mu*( ((r^2 - 3*x^2)/r^5) + ((Ri^2*J2)/2)*(((15*(x^2 + z^2))/r^7) - ((105*x^2*z^2)/r^9) - (3/r^5)) );
delYddDelY = -mu*( ((r^2 - 3*y^2)/r^5) + ((Ri^2*J2)/2)*(((15*(y^2 + z^2))/r^7) - ((105*y^2*z^2)/r^9) - (3/r^5)) );

delXddDelY = mu*( ((3*x*y)/r^5) - ((Ri^2*J2*x*y)/2)*((15/r^7) - ((105*z^2)/r^9)) );
delYddDelX = delXddDelY;

delXddDelZ = mu*( ((3*x*z)/r^5) - ((Ri^2*J2*x*z)/2)*((45/r^7) - ((105*z^2)/r^9)) );
delYddDelZ = mu*( ((3*y*z)/r^5) - ((Ri^2*J2*y*z)/2)*((45/r^7) - ((105*z^2)/r^9)) );

delZddDelX = delXddDelZ;
delZddDelY = delYddDelZ;

delZddDelZ = -mu*( ((r^2 - 3*z^2)/r^5) + ((Ri^2*J2)/2)*(((90*z^2)/r^7) - ((105*z^4)/r^9) - (9/r^5)) );

delXddDelJ2 = (-mu*x/r^3)*((Ri/r)^2)*((15/2)*(z/r)^2 - (3/2));

delYddDelJ2 = (-mu*y/r^3)*((Ri/r)^2)*((15/2)*(z/r)^2 - (3/2));

delZddDelJ2 = (-mu*z/r^3)*((Ri/r)^2)*((15/2)*(z/r)^2 - (9/2));

    % Create matrix blocks for convenience
dynamicsBlock = [
                    delXddDelX, delXddDelY, delXddDelZ;
                    delYddDelX, delYddDelY, delYddDelZ;
                    delZddDelX, delZddDelY, delZddDelZ
                ];

constantsBlock = [
                    delXddDelJ2;
                    delYddDelJ2;
                    delZddDelJ2
                 ];

    % Construct A matrix
A = [
        zeros(3,3), eye(3);
        dynamicsBlock, zeros(3,3)
    ];

end




