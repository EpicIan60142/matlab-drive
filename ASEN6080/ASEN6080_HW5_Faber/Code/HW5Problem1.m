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
filename = "..\Data\HW5Data_J2.mat";

    % Only generate data once
if true
    generateTruthData_MuJ2(earthConst, orbital, x0Perturb, stations, filename, numOrbits, measNoise, dt);
end

load("..\Data\HW5Data_J2.mat");

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

% input("Press Enter to continue, 'Ctrl-C' to exit");

%% Problem 1b. Compare SRIF to LKF without process noise
fprintf("\n1b. Comparing SRIF to LKF\n")

Q0 = 0*diag(sigAccel^2*ones(1,3));

SRIFRun = runSRIF(X0, x0, P0, Q0, uBar, true, earthConst, stations, X_ref, t_ref, 1, true);
LKFRun = runLKF_SNC(X0, x0, P0, Q0, earthConst, stations, X_ref, t_ref, 1, true);

    % Plot difference between SRIF and LKF
titleText = sprintf("SRIF vs. LKF State Difference - no process noise (X_{SRIF} - X_{LKF})");
xLabel = sprintf("Time [sec]"); 
yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
          "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]"];
plotStateError(SRIFRun.t_SRIF, SRIFRun.X_SRIF' - LKFRun.X_LKF', [], [], [], titleText, xLabel, yLabel);

% return

%% Problem 1c. Run SRIF without forcing Rbar to be upper triangular
fprintf("\n1c. Running SRIF without forcing Rbar to be upper triangular\n")

Q0 = 0*diag(sigAccel^2*ones(1,3));

SRIFRunNotTriangular = runSRIF(X0, x0, P0, Q0, uBar, false, earthConst, stations, X_ref, t_ref, 1, true);

    % Plot difference between SRIF and SRIF not triangular
titleText = sprintf("Normal SRIF vs. SRIF w/o upper triangular Rbar");
xLabel = sprintf("Time [sec]"); 
yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
          "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]"];
plotStateError(SRIFRun.t_SRIF, SRIFRun.X_SRIF' - SRIFRunNotTriangular.X_SRIF', [], [], [], titleText, xLabel, yLabel);

