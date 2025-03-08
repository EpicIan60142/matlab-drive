%% ASEN 6080 HW 5 Problem 1 Main Script
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

    % ODE45 options
opt = odeset('AbsTol',1e-12,'RelTol',1e-12);

%% Make truth data
x0Perturb = 0*0.5*[1 1 1 1e-3 1e-3 1e-3]';
numOrbits = 15;
measNoise = true;
dt = 10; % sec
filename = "..\Data\HW5Data_J2J3.mat";

    % Only generate data once
if true
    generateTruthData_MuJ2J3(earthConst, orbital, x0Perturb, stations, filename, numOrbits, measNoise, dt);
end

load("..\Data\HW5Data_J2J3.mat");

titleText = sprintf("Measurement data");
xLabel = "Time [sec]"; 
yLabel = ["\rho [km]", "\rhoDot [km/s]", "Elevation Angle [deg]"];
plotMeasurements(stations, titleText, xLabel, yLabel);

%% Problem 1a. Filter setup
sigR = 1; % km
sigV = 1e-3; % km/s

x0 = 0*0.5*[1 1 1 1e-3 1e-3 1e-3]';
P0 = diag([sigR^2, sigR^2, sigR^2, sigV^2, sigV^2, sigV^2]);

uBar = zeros(3,1);

sigAccel = 1e-8; % Optimal SNC sigma from HW 3

input("Press Enter to continue, 'Ctrl-C' to exit");

%% Problem 1a. Compare SRIF to LKF without process noise

Q0 = 0*diag(sigAccel^2*ones(1,3));

SRIFOut = SRIF(X0, stations, earthConst, P0, x0, Q0, uBar, true);



