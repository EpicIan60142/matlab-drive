%% ASEN 3128 Problem 05 
%   Section 011 - Ian Faber, Sam Mohnacs, Blake Wilson, Kyle Elligott

%% Housekeeping
clc; clear; close all;

%% Constants
g = 9.81; % m/s^2
m = 0.068; % kg
d = 0.06; % m
Km = 0.0024; % N*m/N
Ix = 5.8*10^-5; % kg*m^2
Iy = 7.2*10^-5; % kg*m^2
Iz = 1.0*10^-4; % kg*m^2
nu = 10^-3; % N/(m/2)^2
mu = 2*10^-6; % N*m/(rad/s)^2

%% Problem 5 Root Locus - Longitudinal

% From problem 1
K1_long = 0.001584;
K2_long = 0.00288;
figure
hold on
title("Root Locus of Longitudinal Dynamics")
xlabel("Real")
ylabel("Imaginary")
for K3_long = -0.01:0.0005:0.01
    A = [
            0, -g, 0;
            0, 0, 1;
            -K3_long/Ix, -K2_long/Ix, -K1_long/Ix
        ];

    vals = eig(A);
    
%     k1Vec = [k1Vec; vals(1)]

    x = real(vals);
    y = imag(vals);

    scatter(x,y, 'bo')

end

xline(-0.8)
K3_long = -0.0005;
A = [
            0, -g, 0;
            0, 0, 1;
            -K3_long/Ix, -K2_long/Ix, -K1_long/Ix
    ];

vals = eig(A);
    
%     k1Vec = [k1Vec; vals(1)]

x = real(vals);
y = imag(vals);

scatter(x,y, 'ro')
hold off

%% Problem 5 Root Locus - Lateral

% From problem 1
K1_lat = 0.001276;
K2_lat = 0.00232;
figure
hold on
title("Root Locus of Lateral Dynamics")
xlabel("Real")
ylabel("Imaginary")
for K3_lat = -0.005:0.0005:0.005
    A = [
            0, g, 0;
            0, 0, 1;
            -K3_lat/Ix, -K2_lat/Ix, -K1_lat/Ix
        ];

    vals = eig(A);

    x = real(vals);
    y = imag(vals);

    scatter(x,y, 'bo')

end

xline(-0.8)
K3_lat = 0.0005;
A = [
        0, g, 0;
        0, 0, 1;
        -K3_lat/Ix, -K2_lat/Ix, -K1_lat/Ix
    ];

vals = eig(A);

x = real(vals);
y = imag(vals);

scatter(x,y, 'ro')