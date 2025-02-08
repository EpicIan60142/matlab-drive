function [dX, r, rDot] = calculateOrbit(X, radius)
% Calculates the inertial position and velocity vectors expressed in 
% inertial components for a circular orbit
%
%   Inputs:
%       - X: State vector at a given point in time
%               [Omega; inclination; theta; w_1; w_2; w_3]
%       - radius: Radius of circular orbit
%   Outputs:
%       - dX: Rate of change vector based on the current state
%               [rDot; 
%

EA = X(1:3);
w = X(4:6);

rVec = radius*[1; 0; 0];
rDotVec = radius*w(3)*[0; 1; 0]; 

ON = EA2DCM(EA, [3,1,3]);
NO = ON';

r = NO*rVec;
rDot = NO*rDotVec;

inc = EA(2);
theta = EA(3);

if inc ~= 0
    EADot = (1/sin(inc))*[
                            sin(theta),             cos(theta),             0;
                            cos(theta)*sin(inc),    -sin(theta)*sin(inc),   0;
                            -sin(theta)*cos(inc),   -cos(theta)*cos(inc),   sin(inc);
                         ]*w;
else
    EADot = [0; 0; w(3)];
end
wDot = zeros(3,1);

dX = [EADot; wDot];


end