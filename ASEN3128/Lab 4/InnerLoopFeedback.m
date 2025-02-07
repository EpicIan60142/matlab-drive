function [Fc,Gc] = InnerLoopFeedback(var)
% Section 011 - Ian Faber, Sam Mohnacs, Blake Wilson, Kyle Eligott
%   Function that calculates control forces and moments

m = 0.068; % kg
g = 9.81; % m/s^2

phi = var(4);
theta = var(5);
p = var(10);
q = var(11);

[K1_lat, K2_lat] = CalcGainLat(10);

[K1_long, K2_long] = CalcGainLong(10);

deltaL = -K1_lat*p - K2_lat*phi;
deltaM = -K1_long*q - K2_long*theta;

Gc = [deltaL; deltaM; 0];

Fc = [0; 0; -m*g];

end