%% ASEN 5090 HW 3 Problem 2 Main Script
% By: Ian Faber, 09/21/2025

%% Housekeeping
clc; clear; close all;

addpath("..\utilities\");
addpath("HW3_Code\");
addpath("HW3_Data\");
addpath("..\HW2\HW2_Code\");
addpath("..\HW2\HW2_DATA\");

%% Setup
    % Choose a satellite
PRN = 5;
    
    % Speed of light
c = 299792458; % m/s

%% Part a. Load ephemeris data into an array
    % Read ephemeris data
ephData = read_clean_GPSbroadcast('brdc2230.25n', true);

    % Extract only the satellite of interest
PRNIdx = ephData(:,1) == PRN;

ephData = ephData(PRNIdx, :);

%% Part b. Compute position of GPS satellite based on ephemeris data for the times in the sp3 file from HW 2
    % Read sp3 data
sp3Data = read_sp3("IGS0OPSFIN_20252230000_01D_15M_ORB.SP3");

    % Extract data for satellite of interest
PRNIdx = sp3Data(:,3) == PRN;

sp3Data = sp3Data(PRNIdx, :);

    % Collate times for position calculation
WN = sp3Data(:, 1);
TOW = sp3Data(:, 2);

    % Compute position of GPS satellite from ephemeris data and extract sp3
    % positions
[~, ephPos, ephVel, clkCorr, ~, ~] = eph2pvt2025(ephData, [WN, TOW], PRN);
sp3Pos = sp3Data(:, 4:6)*1000;

TOW = (TOW - TOW(1))/3600;

    % Plot sp3 and ephemeris positions and their difference
figure; tl = tiledlayout(3,1); ax = [];
title(tl, sprintf("Ephemeris and sp3 positions for PRN %.0f", PRN));
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(TOW, ephPos(:,1), 'b-');
    plot(TOW, sp3Pos(:,1), 'r--');
    xlabel("TOW [hr]"); ylabel("X [m]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(TOW, ephPos(:,2), 'b-');
    plot(TOW, sp3Pos(:,2), 'r--');
    xlabel("TOW [hr]"); ylabel("Y [m]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(TOW, ephPos(:,3), 'b-');
    plot(TOW, sp3Pos(:,3), 'r--');
    xlabel("TOW [hr]"); ylabel("Z [m]")
linkaxes(ax, 'x');

figure; tl = tiledlayout(3,1); ax = [];
title(tl, sprintf("Difference between Ephemeris and sp3 Positions for PRN %.0f", PRN))
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(TOW, ephPos(:,1) - sp3Pos(:,1), 'b-');
    xlabel("TOW [hr]"); ylabel("\DeltaX [m]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(TOW, ephPos(:,2) - sp3Pos(:,2), 'b-');
    xlabel("TOW [hr]"); ylabel("\DeltaY [m]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(TOW, ephPos(:,3) - sp3Pos(:,3), 'b-');
    xlabel("TOW [hr]"); ylabel("\DeltaZ [m]")
linkaxes(ax, 'x');

%% Part c. Plot the clock bias based on the a0 and a1 parameters vs. time

figure;
hold on; grid on;
title(sprintf("Satellite clock bias vs. time for PRN %.0f", PRN))
plot(TOW, clkCorr)
xlabel("TOW [hr]"); ylabel("Clock bias [m]")

%% Part d. Calculate and plot relativistic correction vs. time

t_r = 2*dot(ephPos, ephVel, 2)/c^2;
dx_r = c*t_r;

figure;
hold on; grid on;
title(sprintf("Relativistic correction vs. time for PRN %.0f", PRN))
plot(TOW, dx_r, 'b');
xlabel("TOW [hr]"); ylabel("\Deltax_r [m]")



