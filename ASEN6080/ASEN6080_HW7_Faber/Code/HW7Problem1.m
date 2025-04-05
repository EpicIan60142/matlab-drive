%% ASEN 6080 HW 7 Problem 1 Main Script
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
filename = "..\Data\HW6Data_J2J3.mat";

    % Only generate data once
if true
    generateTruthData_MuJ2J3(earthConst, orbital, x0Perturb, stations, filename, numOrbits, measNoise, dt);
end

load("..\Data\HW6Data_J2J3.mat"); % Source of X0, updated stations, X_ref, and t_ref

titleText = sprintf("Measurement data");
xLabel = "Time [sec]"; 
yLabel = ["\rho [km]", "\rhoDot [km/s]", "Elevation Angle [deg]"];
plotMeasurements(stations, titleText, xLabel, yLabel);

input("Press 'Enter' to continue")

%% Problem 1a. Filter setup
    % Covariance matrices
Pxx0 = 1e4*eye(6);
Pcc0 = earthConst.J3^2;

    % Consider parameter setup
C = 0;
c = earthConst.J3 - C;

    % Plotting
plot = [true; true; false; false; false; false; true; true]; % Only plot residuals and state errors

%% Problem 1b/c.




