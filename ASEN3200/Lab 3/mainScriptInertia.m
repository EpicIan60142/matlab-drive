%% ASEN 3200 Lab A-3 Inertia Script
% Group 25 - Ian Faber, John Davis, Steven Liu, Ben Slama

%% Housekeeping
clc; clear; close all;

%% Constants

%% Data Extraction
torqueData = readmatrix("Data\torque10mNm");

time = torqueData(:,1)/1000; % ms to s
cmdTorque = torqueData(:,2)/1000; % mNm to Nm
angV = -torqueData(:,3)*(pi/30); % rpm to rad/s
current = torqueData(:,4); % A

%% Analysis

index = time >= 1 & time <= 6;

[coef, approxCurve] = leastSquares(time(index), angV(index), 1);

alpha = coef(1)

I = mean(cmdTorque(index))/alpha; % kgm^2

%% Plotting
figure
hold on
title("Angular velocity vs. time")
plot(time, angV, 'b.')
plot(time, approxCurve(time), 'k--')
xlabel("Time (sec)")
ylabel("Angular velocity (rad/sec)")

legend("Raw Data", "Best fit")


