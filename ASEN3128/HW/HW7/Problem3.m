%% Housekeeping
clc; clear; close all;

%% Constants
Alat =  [
            -0.0561,     0,          -775.0000,      32.2000
            -0.0038,     -0.4138,    0.4231,         0
            0.0011,      -0.0064,    -0.1456,        0
            0,           1.0000,     0,              0
        ];

Blat = [
            0,          5.601
            -0.14120,   0.1201
            0.00375,    -0.4898
            0,          0
       ];

u0 = -Alat(1,3);

%% Problem 3a

% C = zeros(size(Alat));
% C(1,1) = 1/u0;
C = [1/u0, 0, 0, 0];
D = 0;

[num, den] = ss2tf(Alat, Blat(:,2), C, D);

tf = tf(num, den);

%% Problem 3b
del_r = deg2rad(2); % 2 degree rudder step input 

[y, t] = step(tf);

y = del_r*y;

figure
title("Sideslip step response to rudder input")
plot(t, rad2deg(y));
xlabel("Time (sec)")
ylabel("Sideslip (deg)")

%% Problem 3c
eigVals = roots(den);

dutch = eigVals(1:2);
real = real(dutch(1));
imag = imag(dutch(1));

wn = sqrt(real^2 + imag^2);
zeta = -real/wn;

