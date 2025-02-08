%% ASEN 3112 Lab 1 Problem 1 Script (Closed Thin Wall Analysis)
% Section 013 Group 6 - Liam Casey, Ian Faber, Chayce Johnson, Tori Morgheim, Ethan Stapp
% Code by: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Constants
dExt = 0.75; % Exterior diameter, inches
rExt = dExt/2; % Exterior radius, inches
t = 1/16; % Wall thickness, inches
L = 1; % Extensometer gauge length, inches
G = 3.75 * 10^6; % Shear modulus, psi

%% Data Extraction
CTWData = readmatrix('Data\CTW400inlbf_05.txt');

time = CTWData(:,1);
torsionalAngle = CTWData(:,2);
shearStrain = CTWData(:,3);
torque = CTWData(:,4);
axialStrain = CTWData(:,5);

%% Calculate shear strain from torsional Angle
% phi = gamma*L/Re
gamma = torsionalAngle*rExt/L;

%% Calculate Least Squares fit of both shear strain measurements
[coef1, approxCurve1] = leastSquares(torque, shearStrain, 1);
GJ_1 = coef1(1);
label1 = sprintf("Line of best fit: \n \\gamma = %.3e*T + %.3e", coef1(1), coef1(2));

[coef2, approxCurve2] = leastSquares(torque, gamma, 1);
GJ_2 = coef2(1);
label2 = sprintf("Line of best fit: \n \\gamma = %.3e*T + %.3e", coef2(1), coef2(2));

%% Plotting
figure
hold on;
title("CTW Extensometer Shear Strain vs. Torque")
plot(torque, shearStrain, 'b.', 'MarkerSize', 1);
plot(torque, approxCurve1(torque), 'k--')
xlabel("Torque (in-lbf)")
ylabel("Shear Strain (deg)")
legend("Raw extensometer data", label1, 'Location', 'best')

figure
hold on;
title("CTW Calculated Shear Strain vs. Torque")
plot(torque, gamma, 'r.')
plot(torque, approxCurve2(torque), 'k--')
xlabel("Torque (in-lbf)")
ylabel("Calculated Shear Strain (deg)")
legend("Calculated shear strain", label2, 'Location', 'best')

%% Run problem 2 script if desired
run('Problem2.m')

