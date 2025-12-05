%% ASEN 5254 HW 9 Main
% By: Ian Faber, 10/31/2025

%% Housekeeping
clc; clear; close all;

%% Plot control data
data = readmatrix("controlData.csv");

time = data(:,1);
velocity = data(:,2);
steeringAngle = data(:,3);

figure(); tl = tiledlayout(2,1); ax = [];
title(tl, "Parallel Parking Control Inputs")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    title("Velocity vs. time")
    plot(time, velocity, 'b');
    xlabel("Time [sec]"); ylabel("Velocity")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    title("Steering Angle vs. time")
    plot(time, steeringAngle, 'r');
    xlabel("Time [sec]"); ylabel("Steering Angle [rad]")
linkaxes(ax, 'x');