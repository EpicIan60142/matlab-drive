%% ASEN 5090 HW 2 Problem 1 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

addpath("HW2_Code\")
addpath("HW2_DATA\")

%% Setup
sp3Data = read_sp3("IGS0OPSFIN_20252230000_01D_15M_ORB.SP3");
almanacData = read_GPSyuma("YUMA223.alm", 2);
PRN = 16; % Satellite to pull data for

%% Part a/b. Plot SP3 and almanac data vs. time
    % Choose rows for this PRN
idxSat_SP3 = sp3Data(:,3) == PRN;

    % Pull out SP3 data
t_SP3 = sp3Data(idxSat_SP3,2);
t_SP3 = (t_SP3 - t_SP3(1))/3600; % Subtract initial time and convert to hours
X_SP3 = sp3Data(idxSat_SP3,4);
Y_SP3 = sp3Data(idxSat_SP3,5);
Z_SP3 = sp3Data(idxSat_SP3,6);
SP3Pos = [X_SP3, Y_SP3, Z_SP3]*1000; % Convert from km to m

WN = sp3Data(idxSat_SP3, 1);
TOW = sp3Data(idxSat_SP3, 2);

    % Pull out almanac data
[~, almanacPos] = alm2pos(almanacData, [WN, TOW], PRN);

    % Plot data
yLabels = ["X [m]", "Y [m]", "Z [m]"];
fig = figure(); tl = tiledlayout(3,1); ax = [];
title(tl, sprintf("GPS Satellite ECEF positions for PRN %.0f", PRN))
for k = 1:3
    nt = nexttile;
        hold on; grid on;
        sp3Plot = plot(t_SP3, SP3Pos(:,k), 'b-');
        almPlot = plot(t_SP3, almanacPos(:,k), 'r--');
        xlabel("TOW [hr]"); ylabel(yLabels(k))
    ax = [ax; nt];
end
linkaxes(ax, 'x')
legend([sp3Plot, almPlot], ["SP3 Data", "YUMA Almanac Data"], 'location', 'northeastoutside');

%% Part c. Plot the differences between the YUMA Almanac and SP3 Positions
    % Calculate difference
difference = almanacPos - SP3Pos;

    % Plot difference
yLabels = ["\DeltaX [m]", "\DeltaY [m]", "\DeltaZ [m]"];
fig = figure(); tl = tiledlayout(3,1); ax = [];
title(tl, sprintf("YUMA Almanac minus SP3 ECEF Positions for PRN %.0f", PRN))
for k = 1:3
    nt = nexttile;
        hold on; grid on;
        plot(t_SP3, difference(:,k), 'b-');
        xlabel("TOW [hr]"); ylabel(yLabels(k))
    ax = [ax; nt];
end
linkaxes(ax, 'x')




