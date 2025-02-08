%% ASEN 6080 HW 2 Problem 2 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath(genpath("..\..\Utilities\"));

    % Extract parameters for Earth
earthConst = getPlanetConst();

    % Set up stations
stations = makeStations(earthConst);

    % Extract orbital elements
orbital = getOrbitConst();

%% Problem 2a: Make Truth Data with J3
x0Perturb = 1*0.5*[1 1 1 1e-3 1e-3 1e-3]';
numOrbits = 15;
dt = 10; % sec
filename = "HW2Problem2Data.mat";

    % Only generate data once
if false
    generateTruthData_MuJ2J3(earthConst, orbital, x0Perturb, stations, filename, numOrbits, dt);
end

J3Data = load("HW2Problem2Data.mat");

    % Generate J2 data from HW 1
J2Data = generateTruthData_MuJ2(earthConst, orbital, x0Perturb, stations, [], numOrbits, dt);

    % Plot differences in truth data
titleText = sprintf("State differences between J3 and J2 data");
xLabel = "Time [sec]"; 
yLabel = ["\DeltaX [km]", "\DeltaY [km]", "\DeltaZ [km]", ...
          "\DeltaXdot [km/s]", "\DeltaYdot [km/s]", "\DeltaZdot [km/s]"];
plotStateError(J3Data.t_ref, J3Data.X_ref - J2Data.X_ref, [], [], titleText, xLabel, yLabel);

titleText = sprintf("Measurement differences between J3 and J2 data");
xLabel = "Time [sec]"; 
yLabel = ["\Delta\rho [km]", "\Delta\rhoDot [km/s]", "\DeltaElevation Angle [deg]"];
for k = 1:length(J3Data.stations)
    stationComp(k).tMeas = J3Data.stations(k).tMeas;
    stationComp(k).rho = J3Data.stations(k).rho - J2Data.stations(k).rho;
    stationComp(k).rhoDot = J3Data.stations(k).rhoDot - J2Data.stations(k).rhoDot;
    stationComp(k).elAngle = J3Data.stations(k).elAngle - J2Data.stations(k).elAngle;
    stationComp(k).color = J3Data.stations(k).color;
end
plotMeasurements(stationComp, titleText, xLabel, yLabel);

%% Problem 2b: Filter setup
sigR = 1; % km
sigV = 1e-3; % km/s

x0 = 1*0.5*[1 1 1 1e-3 1e-3 1e-3]';
P0 = diag([sigR^2, sigR^2, sigR^2, sigV^2, sigV^2, sigV^2]);

%% Problem 2b: Run filters on new data
    
    % Batch
batchRun = runBatch(J3Data.X0, x0, P0, earthConst, J3Data.stations, J3Data.X_ref, J3Data.tspan, J3Data.opt);

    % LKF
LKFRun = runLKF(J3Data.X0, x0, P0, earthConst, J3Data.stations, J3Data.X_ref, J3Data.t_ref);

    % EKF
numMeas = 10;
X0_EKF = LKFRun.X_LKF(:,numMeas);
P0_EKF = LKFRun.LKFOut.PEst{numMeas};
t_EKF_start = LKFRun.t_LKF(numMeas);
EKFRun = runEKF(X0_EKF, P0_EKF, earthConst, J3Data.stations, J3Data.X_ref, J3Data.t_ref, t_EKF_start);

