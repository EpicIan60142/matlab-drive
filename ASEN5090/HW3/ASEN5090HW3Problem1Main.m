%% ASEN 5090 HW 3 Problem 1 Main Script
% By: Ian Faber, 09/21/2025

%% Housekeeping
clc; clear; close all;

addpath("..\utilities\")
addpath("HW3_Code\")
addpath("HW3_Data\")

%% Part a. Read RINEX file and choose a PRN with a long pass over Boulder
    % Read data with built-in MATLAB function
rinexData = rinexread("NIST00USA_R_20252230000_01D_30S_MO.rnx");

GPSData = rinexData.GPS;

    % Choose a satellite PRN
PRN = 5;

    % Extract only data for that satellite
PRNIdx = GPSData.SatelliteID == PRN;

GPSData = GPSData(PRNIdx, :);

%% Part b. Plot C1C pseudorange vs. time for the chosen PRN
figure;
hold on; grid on;
title(sprintf("C1C Pseudorange vs. time for PRN %.0f", PRN));
plot(GPSData.Time, GPSData.C1C, 'b.');
xlabel("Time"); ylabel("C1C Pseudorange [m]");

%% Part c. Plot S1C signal to noise ratio vs. time for the chosen PRN
figure;
hold on; grid on;
title(sprintf("S1C SNR vs. time for PRN %.0f", PRN));
plot(GPSData.Time, GPSData.S1C, 'b.');
xlabel("Time"); ylabel("S1C SNR [dB-Hz]");

%% Part d. Choose another pseudorange or carrier phase measurement and plot against time
    % GPS wavelength and speed of light
c = 299792458; % m/s
lambda = c/(1575.42e6); % m

L1C_m = lambda*GPSData.L1C;

figure;
hold on; grid on;
title(sprintf("L1C carrier phase (in meters) vs. time for PRN %.0f", PRN));
plot(GPSData.Time, L1C_m, 'b.');
xlabel("Time"); ylabel("L1C carrier phase [m]");

figure; tl = tiledlayout(1,3);
title(tl, sprintf("Difference in C1C and L1C pseudorange for PRN %.0f", PRN))
nexttile([1 2]);
    hold on; grid on;
    title("C1C - L1C vs. time")
    plot(GPSData.Time, GPSData.C1C - L1C_m, 'b.');
    xlabel("Time"); ylabel("\Delta\rho [m]");
    yLim = ylim();
nexttile;
    hold on; grid on;
    title("Range Difference Histogram")
    histogram(GPSData.C1C - L1C_m, 'FaceColor', 'b', 'Orientation', 'horizontal')
    xlabel("Bin Count"); ylim(yLim);






