function Htilde = MeasurementPartials_RngRngRate_stat(X,X_s)
% Function that outputs the measurement partials matrix for orbital
% measurements using range and range rate for ground stations
%   Inputs:
%       - X: Spacecraft state arranged as follows:
%            [x; y; z; xDot; yDot; zDot; J2]
%       - X_s: Station state arranged as follows:
%            [x_s; y_s; z_s; xDot_s; yDot_s; zDot_s]
%
%   Outputs:
%       - Htilde: Measurement partials matrix
%
%   By: Ian Faber, 01/24/2025
%

x = X(1);
y = X(2);
z = X(3);
xDot = X(4);
yDot = X(5);
zDot = X(6);

x_s = X_s(1);
y_s = X_s(2);
z_s = X_s(3);
xDot_s = X_s(4);
yDot_s = X_s(5);
zDot_s = X_s(6);

rho = norm(X(1:3) - X_s(1:3));

delRhoDelXs =  (x_s-x)/rho;
delRhoDelYs =  (y_s-y)/rho;
delRhoDelZs =  (z_s-z)/rho;

delRhoDotDelXs = (rho^2*(xDot_s - xDot) + (x-x_s)*((x-x_s)*(xDot-xDot_s) + (y-y_s)*(yDot-yDot_s) + (z-z_s)*(zDot-zDot_s)))/(rho^3);
delRhoDotDelYs = (rho^2*(yDot_s - yDot) + (y-y_s)*((x-x_s)*(xDot-xDot_s) + (y-y_s)*(yDot-yDot_s) + (z-z_s)*(zDot-zDot_s)))/(rho^3);
delRhoDotDelZs = (rho^2*(zDot_s - zDot) + (z-z_s)*((x-x_s)*(xDot-xDot_s) + (y-y_s)*(yDot-yDot_s) + (z-z_s)*(zDot-zDot_s)))/(rho^3);
 
Htilde = [
            delRhoDelXs, delRhoDelYs, delRhoDelZs;
            delRhoDotDelXs, delRhoDotDelYs, delRhoDotDelZs
         ];

% R = X(1:3);
% Rs = X_s(1:3);
% V = X(4:6);
% 
% wTilde = tilde([0;0;7.2921158553e-5]);
% rhoVec = R - Rs;
% rho = norm(rhoVec);
% 
% stationPosPart = -rhoVec'/rho;
% stationVelPart = ((rho^2)*(wTilde*R - V) + rhoVec'*(V - wTilde*Rs)*rhoVec)/rho^3;
% 
% Htilde = [stationPosPart;stationVelPart'];

end