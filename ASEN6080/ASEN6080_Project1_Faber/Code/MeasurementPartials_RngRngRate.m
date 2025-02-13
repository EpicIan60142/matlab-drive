function Htilde = MeasurementPartials_RngRngRate(X, statVis, pConst)
% Function that outputs the measurement partials matrix for orbital
% measurements using range and range rate for spacecraft
%   Inputs:
%       - X: Problem state vector in km:
%            X = [x; y; z; xDot; yDot; zDot; mu; J2; Cd; stationPos_1;
%                 stationPos_2; stationPos_3]
%       - statVis: Which station index produced this measurement
%       - pConst: Planetary constants vector as defined by getPlanetConst.m
%
%   Outputs:
%       - Htilde: Measurement partials matrix
%
%   By: Ian Faber, 02/11/2025
%
    % Spacecraft state block

    % Extract spacecraft state
x = X(1);
y = X(2);
z = X(3);
xDot = X(4);
yDot = X(5);
zDot = X(6);

    % Extract visible station state
idx = statVis;
x_s = X(7 + 3*idx);
y_s = X(7 + 3*idx + 1);
z_s = X(7 + 3*idx + 2);

    % Calculate visible station velocity
v_s = cross(pConst.wPlanet, [x_s; y_s; z_s]);
xDot_s = v_s(1);
yDot_s = v_s(2);
zDot_s = v_s(3);

R = [x; y; z];
Rs = [x_s; y_s; z_s];
rho = norm(R - Rs);

delRhoDelX =  (x-x_s)/rho;
delRhoDelY =  (y-y_s)/rho;
delRhoDelZ =  (z-z_s)/rho;

delRhoDotDelX = (rho^2*(xDot - xDot_s) - (x-x_s)*((x-x_s)*(xDot-xDot_s) + (y-y_s)*(yDot-yDot_s) + (z-z_s)*(zDot-zDot_s)))/(rho^3);
delRhoDotDelY = (rho^2*(yDot - yDot_s) - (y-y_s)*((x-x_s)*(xDot-xDot_s) + (y-y_s)*(yDot-yDot_s) + (z-z_s)*(zDot-zDot_s)))/(rho^3);
delRhoDotDelZ = (rho^2*(zDot - zDot_s) - (z-z_s)*((x-x_s)*(xDot-xDot_s) + (y-y_s)*(yDot-yDot_s) + (z-z_s)*(zDot-zDot_s)))/(rho^3);

delRhoDotDelXDot = delRhoDelX;
delRhoDotDelYDot = delRhoDelY;
delRhoDotDelZDot = delRhoDelZ;

HBlock_1 = [
                    delRhoDelX, delRhoDelY, delRhoDelZ, zeros(1,6);
                    delRhoDotDelX, delRhoDotDelY, delRhoDotDelZ, delRhoDotDelXDot, delRhoDotDelYDot, delRhoDotDelZDot, zeros(1,3)
                 ];

    % Station state block
R = [x; y; z];
V = [xDot; yDot; zDot];
wTilde = tilde(pConst.wPlanet);

HBlock_2 = [];
for k = 1:3
    Rs = [
            X(7 + 3*k);
            X(7 + 3*k + 1);
            X(7 + 3*k + 2);
         ];

    rhoVec = R - Rs;
    rho = norm(rhoVec);

    stationPosPart = -rhoVec'/rho;
    stationVelPart = ((rho^2)*(wTilde*R - V) + rhoVec'*(V - wTilde*Rs)*rhoVec)/rho^3;

    HBlock_2 = [HBlock_2, [stationPosPart; stationVelPart']];
end

Htilde = [HBlock_1, HBlock_2];

end