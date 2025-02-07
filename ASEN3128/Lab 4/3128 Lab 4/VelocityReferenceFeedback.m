function [Fc, Gc] = VelocityReferenceFeedback(t,var)
% Section 011 - Ian Faber, Sam Mohnacs, Blake Wilson, Kyle Eligott
%   Function that calculates control forces and moments based on velocity
%   reference

m = 0.068; % kg
g = 9.81; % m/s^2

t1 = 0;
t2 = 2.005;

wr = 0;

if t > t1 && t < t2
    vr = 1*0.5; % Multiply by 1 to test lateral controller
    ur = 0*0.5; % Multiply by 1 to test longitudinal controller
else
    vr = 0;
    ur = 0;
end

phi = var(4);
theta = var(5);
u = var(7);
v = var(8);
w = var(9);
p = var(10);
q = var(11);
r = var(12);

K1_lat = 0.001276;
K2_lat = 0.00232;
K3_lat = 0.0005;

K1_long = 0.001584;
K2_long = 0.00288;
K3_long = -0.0005;

deltaL = -K1_lat*p - K2_lat*phi + K3_lat*(vr - v);
deltaM = -K1_long*q - K2_long*theta + K3_long*(ur - u);

deltaN = -0.004*r;

Gc = [deltaL; deltaM; deltaN];

Fc = [0; 0; -m*g + 0*0.4*(wr - w)]; % Gain on control force for fun, NOT 
                                    % part of assignment (turn on by 
                                    % multiplying by 1)

end