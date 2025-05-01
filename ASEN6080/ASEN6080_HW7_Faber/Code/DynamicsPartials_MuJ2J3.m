function A = DynamicsPartials_MuJ2J3(X, mu, J2, J3, Ri)
% Function that outputs the acceleration partials matrix for orbital 
% dynamics including contributions from Mu, J2, and J3
%   Inputs:
%       - X: System state in SI units 
%           [(x,y,z) -> km, (xDot, yDot, zDot) -> km/s
%           X = [x; y; z; xDot; yDot; zDot; mu; J2; J3]
%       - Ri: Reference radius [km]
%
%   Outputs:
%       - A: Acceleration partials matrix
%
%   By: Ian Faber, 04/08/2025
%

    % Extract states
x = X(1);
y = X(2);
z = X(3);

    % Define range
r = sqrt(x^2 + y^2 + z^2);

    % Define non-trivial partials
delXddDelX = -mu*( ((r^2 - 3*x^2)/r^5) + ((Ri^2*J2)/2)*(((15*(x^2 + z^2))/r^7) - ((105*x^2*z^2)/r^9) - (3/r^5)) + ((Ri^3*J3*z)/2)*(((35*(z^2 + (3*x^2)))/r^9) - ((315*x^2*z^2)/r^11) - (15/r^7)) );
delYddDelY = -mu*( ((r^2 - 3*y^2)/r^5) + ((Ri^2*J2)/2)*(((15*(y^2 + z^2))/r^7) - ((105*y^2*z^2)/r^9) - (3/r^5)) + ((Ri^3*J3*z)/2)*(((35*(z^2 + (3*y^2)))/r^9) - ((315*y^2*z^2)/r^11) - (15/r^7)) );

delXddDelY = mu*( ((3*x*y)/r^5) - ((Ri^2*J2*x*y)/2)*((15/r^7) - ((105*z^2)/r^9)) - ((Ri^3*J3*x*y)/2)*(((105*z)/r^9) - ((315*z^3)/r^11)) );
delYddDelX = delXddDelY;

delXddDelZ = mu*( ((3*x*z)/r^5) - ((Ri^2*J2*x*z)/2)*((45/r^7) - ((105*z^2)/r^9)) - ((Ri^3*J3*x)/2)*(((210*z^2)/r^9) - ((315*z^4)/r^11) - (15/r^7)) );
delYddDelZ = mu*( ((3*y*z)/r^5) - ((Ri^2*J2*y*z)/2)*((45/r^7) - ((105*z^2)/r^9)) - ((Ri^3*J3*y)/2)*(((210*z^2)/r^9) - ((315*z^4)/r^11) - (15/r^7)) );

delZddDelX = delXddDelZ;
delZddDelY = delYddDelZ;

delZddDelZ = -mu*( ((r^2 - 3*z^2)/r^5) + ((Ri^2*J2)/2)*(((90*z^2)/r^7) - ((105*z^4)/r^9) - (9/r^5)) + ((Ri^3*J3*z)/2)*(((350*z^2)/r^9) - ((315*z^4)/r^11) - (75/r^7)) );

delXddDelMu = (-x/r^3)*(1 + ((Ri/r)^2)*J2*((15/2)*(z/r)^2 - (3/2)) + ((Ri/r)^3)*J3*((35/2)*(z/r)^3 - (15/2)*(z/r)));
delXddDelJ2 = (-mu*x/r^3)*((Ri/r)^2)*((15/2)*(z/r)^2 - (3/2));
delXddDelJ3 = (-mu*x/r^3)*((Ri/r)^3)*((35/2)*(z/r)^3 - (15/2)*(z/r));

delYddDelMu = (-y/r^3)*(1 + ((Ri/r)^2)*J2*((15/2)*(z/r)^2 - (3/2)) + ((Ri/r)^3)*J3*((35/2)*(z/r)^3 - (15/2)*(z/r)));
delYddDelJ2 = (-mu*y/r^3)*((Ri/r)^2)*((15/2)*(z/r)^2 - (3/2));
delYddDelJ3 = (-mu*y/r^3)*((Ri/r)^3)*((35/2)*(z/r)^3 - (15/2)*(z/r));

delZddDelMu = (-z/r^3)*(1 + ((Ri/r)^2)*J2*((15/2)*(z/r)^2-(9/2))) - (1/r^2)*((Ri/r)^3)*J3*((35/2)*(z/r)^4 - 15*(z/r)^2 + (3/2));
delZddDelJ2 = (-mu*z/r^3)*((Ri/r)^2)*((15/2)*(z/r)^2-(9/2));
delZddDelJ3 = -(mu/r^2)*((Ri/r)^3)*((35/2)*(z/r)^4 - 15*(z/r)^2 + (3/2));

    % Create matrix blocks for convenience
dynamicsBlock = [
                    delXddDelX, delXddDelY, delXddDelZ;
                    delYddDelX, delYddDelY, delYddDelZ;
                    delZddDelX, delZddDelY, delZddDelZ
                ];

constantsBlock = [
                    delXddDelMu, delXddDelJ2, delXddDelJ3;
                    delYddDelMu, delYddDelJ2, delYddDelJ3;
                    delZddDelMu, delZddDelJ2, delZddDelJ3
                 ];

    % Construct A matrix
A = [
        zeros(3,3), eye(3);
        dynamicsBlock, zeros(3,3)
    ];

end




