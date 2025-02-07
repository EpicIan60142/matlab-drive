function [Fc,Gc] = InnerLoopFeedback(var)
% Section 011 - Ian Faber, Sam Mohnacs, Blake Wilson, Kyle Eligott
%   Function that calculates control forces and moments

m = 0.068; % kg
g = 9.81; % m/s^2
n_lat = 10;
n_long = 10;

phi = var(4);
theta = var(5);
p = var(10);
q = var(11);
r = var(12);

[K1_lat, K2_lat] = CalcGainLat(n_lat);

[K1_long, K2_long] = CalcGainLong(n_long);

K1_lat = 0.001276;
K2_lat = 0.00232;

K1_long = 0.001584;
K2_long = 0.00288;

deltaL = -K1_lat*p - K2_lat*phi;
deltaM = -K1_long*q - K2_long*theta;

Gc = [deltaL; deltaM; -0.004*r];

Fc = [0; 0; -m*g];

end