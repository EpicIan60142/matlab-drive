%% Altimeter Data Analysis
% Diabetics of Space Z, code by Ian Faber
% Created: 4/20/22
% Analyzes altimeter data for height and flight time

%% Housekeeping
clc; clear; close all;

%% Constants
rho = 1.14; % kg/m^3
g = 9.81; % m/s^2

%% Data extraction
data = load('JoshuaCampB10_raw_processed_smoothed');

time = data(:,1);
pressureRaw = data(:,2);
temperature = data(:,3);

%% Processing
pressurePa = 100*(pressureRaw - pressureRaw(1));
height = -pressurePa/(g*rho);

%% Plotting
figure()
hold on;
title("Bottle Rocket Height")
plot(time, height);
xlabel("Time (sec)")
ylabel("Height (m)")
hold off;


