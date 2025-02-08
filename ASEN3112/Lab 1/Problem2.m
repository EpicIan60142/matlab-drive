%% ASEN 3112 Lab 1 Problem 2 Script (Open Thin Wall Analysis)
% Section 013 Group 6 - Liam Casey, Ian Faber, Chayce Johnson, Tori Morgheim, Ethan Stapp
% Code by: Ethan Stapp, Ian Faber

%% Housekeeping
clc; clear; close all;

%% Constants
dExt = 0.75; % Exterior diameter, inches
rExt = dExt/2; % Exterior radius, inches
t = 1/16; % Wall thickness, inches
rInt = rExt-t; % Interior radius, inches
L = 1; % Extensometer gauge length, inches
G = 3.75 * 10^6; % Shear modulus, psi
gripLength = 10; % Length between grips of encoder, inches

%% Data Extraction
OTWData = readmatrix('OTW20inlbf_004.txt');

time = OTWData(:,1); % seconds
torsionalAngle = deg2rad(OTWData(:,2)); % Torsional angle, radians
torsionalAngle = torsionalAngle - min(torsionalAngle);
shearStrain = deg2rad(OTWData(:,3)); % Shear strain, rad
shearStrain = shearStrain - min(shearStrain);
torque = OTWData(:,4); % Torque, in-lbs
axialStrain = OTWData(:,5); % Axial strain
dphi_dx_extensometer = torsionalAngle./L; % Dividing the torsional angle by L

%% Calculate shear strain from torsional Angle
% phi = gamma*L/t
dphi_dx_encoder = torsionalAngle./gripLength; % Dividing the torsional angle by the distance between encoder grips
gamma = torsionalAngle * t/L;
gamma = gamma - min(gamma);

%% Calculate Least Squares fit of both shear strain measurements
[coef1, approxCurve1] = leastSquares(dphi_dx_encoder, torque, 1);
GJ_encoder = coef1(1);

[coef2, approxCurve2] = leastSquares(dphi_dx_extensometer, torque, 1);
GJ_extensometer = coef2(1);

% Calculating torsional rigidity using thin wall theory:
b = 2*pi*(rExt+rInt)/2;
J_beta = (1/3) * b * (t^3);
GJ_analytical = G * J_beta;

% Plotting
figure
hold on;
title("OTW Extensometer Shear Strain vs. Torque")
plot(shearStrain, torque, 'b.', 'MarkerSize', 1);
ylabel("Torque (in-lbf)")
xlabel("Gamma (rad)")
legend("Raw extensometer data", 'Location', 'best')

figure
hold on;
title("OTW Calculated Shear Strain vs. Torque")
plot(gamma, torque, 'r.')
ylabel("Torque (in-lbf)")
xlabel("Calculated Gamma (rad)")
legend("Calculated shear strain", 'Location', 'best')

% figure
% hold on;
% plot(axialStrain*L + L, torsionalAngle - min(torsionalAngle), 'g.');
% xlabel("X (in)");
% ylabel("Torsional Angle (rad)")
% title("Torsional Angle vs. X")

%% Calculating Percent Error:
percentError1 = (abs((GJ_extensometer-GJ_analytical))/GJ_analytical).*100;
percentError2 = (abs((GJ_encoder-GJ_analytical))/GJ_encoder).*100;



