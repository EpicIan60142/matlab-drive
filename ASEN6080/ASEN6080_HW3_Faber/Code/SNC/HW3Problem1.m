%% ASEN 6080 HW 3 Problem 1 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath("..\")
addpath(genpath("..\..\..\Utilities"));

    % Extract parameters for Earth
earthConst = getPlanetConst();

    % Set up stations
stations = makeStations(earthConst);

    % Extract orbital elements
orbital = getOrbitConst();

%% Make Truth Data
x0Perturb = 1*0.5*[1 1 1 1e-3 1e-3 1e-3]';
numOrbits = 15;
dt = 10; % sec
filename = "HW3Problem1Data.mat";

    % Only generate data once
if true
    generateTruthData_MuJ2J3(earthConst, orbital, 0*x0Perturb, stations, filename, numOrbits, dt);
end

load("HW3Problem1Data.mat");

titleText = sprintf("Measurement data");
xLabel = "Time [sec]"; 
yLabel = ["\rho [km]", "\rhoDot [km/s]", "Elevation Angle [deg]"];
plotMeasurements(stations, titleText, xLabel, yLabel);

%% Problem 1b: Filter setup
sigR = 1; % km
sigV = 1e-3; % km/s

x0 = 1*0.5*[1 1 1 1e-3 1e-3 1e-3]';
P0 = diag([sigR^2, sigR^2, sigR^2, sigV^2, sigV^2, sigV^2]);

input("Press Enter to continue, 'Ctrl-C' to exit");



