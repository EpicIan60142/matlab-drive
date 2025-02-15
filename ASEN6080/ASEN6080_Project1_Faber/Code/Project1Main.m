%% ASEN 6080 Project 1 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath(genpath("..\..\Utilities"))

earthConst = getPlanetConst();
scConst = getSCConst();
stations = makeStations();

%% Load and parse provided data
data = readmatrix("..\Data\project.txt");

stations = readData(stations, data);

%% Plot measurements
titleText = "Range and Range-Rate Measurements";
xLabel = "Time [sec]"; yLabel = ["Range [m]","Range-Rate [m/s]"];
plotMeasurements(stations, titleText, xLabel, yLabel);

%% Filter setup
stationX0 = [];
for k = 1:length(stations)
    stationX0 = [stationX0; stations(k).X0];
end

X0 = [scConst.X0_cart; earthConst.mu; earthConst.J2; scConst.Cd; stationX0];
P0 = diag([
            1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e20, 1e6, 1e6, ...
            1e-10, 1e-10, 1e-10, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6
          ]);

opt = odeset('AbsTol',1e-12,'RelTol',1e-12);

%% Run Batch - all measurements, original P0
batchRun = runBatch(X0, zeros(size(X0)), P0, earthConst, scConst, stations, data(:,1), [], opt, 10);

%% Run LKF - original P0
LKFRun = runLKF(X0, zeros(size(X0)), P0, earthConst, scConst, stations, 5);


