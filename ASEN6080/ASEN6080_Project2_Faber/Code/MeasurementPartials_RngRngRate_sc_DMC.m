function Htilde = MeasurementPartials_RngRngRate_sc_DMC(X,X_s)
% Function that outputs the measurement partials matrix for orbital
% measurements using range and range rate for spacecraft
%   Inputs:
%       - X: Spacecraft state arranged as follows:
%            [x; y; z; xDot; yDot; zDot; C_R]
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

delRhoDelX =  (x-x_s)/rho;
delRhoDelY =  (y-y_s)/rho;
delRhoDelZ =  (z-z_s)/rho;

delRhoDotDelX = (rho^2*(xDot - xDot_s) - (x-x_s)*((x-x_s)*(xDot-xDot_s) + (y-y_s)*(yDot-yDot_s) + (z-z_s)*(zDot-zDot_s)))/(rho^3);
delRhoDotDelY = (rho^2*(yDot - yDot_s) - (y-y_s)*((x-x_s)*(xDot-xDot_s) + (y-y_s)*(yDot-yDot_s) + (z-z_s)*(zDot-zDot_s)))/(rho^3);
delRhoDotDelZ = (rho^2*(zDot - zDot_s) - (z-z_s)*((x-x_s)*(xDot-xDot_s) + (y-y_s)*(yDot-yDot_s) + (z-z_s)*(zDot-zDot_s)))/(rho^3);

delRhoDotDelXDot = delRhoDelX;
delRhoDotDelYDot = delRhoDelY;
delRhoDotDelZDot = delRhoDelZ;
 
Htilde = [
            delRhoDelX, delRhoDelY, delRhoDelZ, zeros(1,4);
            delRhoDotDelX, delRhoDotDelY, delRhoDotDelZ, delRhoDotDelXDot, delRhoDotDelYDot, delRhoDotDelZDot, 0
         ];

Htilde = [Htilde, zeros(size(Htilde,1),3)];

end