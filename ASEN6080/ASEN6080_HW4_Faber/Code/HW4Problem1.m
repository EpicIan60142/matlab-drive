%% ASEN 6080 HW 3 Problem 1 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath("..\")
addpath(genpath("..\..\..\Utilities\"));

    % Extract parameters for Earth
earthConst = getPlanetConst();

    % Set up stations
stations = makeStations(earthConst);

    % Extract orbital elements
orbital = getOrbitConst();

%% Make Truth Data
x0Perturb = 0*0.5*[1 1 1 1e-3 1e-3 1e-3]';
numOrbits = 15;
measNoise = true;
dt = 10; % sec
filename = "..\Data\HW4Data_J2J3.mat";

    % Only generate data once
if false
    generateTruthData_MuJ2J3(earthConst, orbital, x0Perturb, stations, filename, numOrbits, measNoise, dt);
end

load("..\Data\HW4Data_J2J3.mat");

titleText = sprintf("Measurement data");
xLabel = "Time [sec]"; 
yLabel = ["\rho [km]", "\rhoDot [km/s]", "Elevation Angle [deg]"];
plotMeasurements(stations, titleText, xLabel, yLabel);

%% Problem 1b: Filter setup
sigR = 1; % km
sigV = 1e-3; % km/s

x0 = 0*0.5*[1 1 1 1e-3 1e-3 1e-3]';
P0 = diag([sigR^2, sigR^2, sigR^2, sigV^2, sigV^2, sigV^2]);

sigAccel = 1e-8; % Optimal SNC sigma from HW 3

input("Press Enter to continue, 'Ctrl-C' to exit");


%% Problem 1b: Choose optimal sigma and show state errors
    % Based on plots, sigma = 1e-8 balances both postfit and 3D RMS
fprintf("\nPlotting state errors vs. time for sigma = %.3e km/s^2\n", sigAccel)

    % Define optimal Q0
if measNoise
    Q0 = 1*diag(sigAccel^2*ones(1,3));
else
    Q0 = zeros(3,3);
end

    % Run filter
LKFOpt = runLKF_SNC(X0, x0, P0, Q0, earthConst, stations, X_ref, t_ref, 1, true);






