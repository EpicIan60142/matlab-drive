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

    % ODE45 options
opt = odeset('AbsTol',1e-12,'RelTol',1e-12);

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

%% Problem 1a: Filter setup
sigR = 1; % km
sigV = 1e-3; % km/s

x0 = 0*0.5*[1 1 1 1e-3 1e-3 1e-3]';
P0 = diag([sigR^2, sigR^2, sigR^2, sigV^2, sigV^2, sigV^2]);

sigAccel = 1e-8; % Optimal SNC sigma from HW 3

input("Press Enter to continue, 'Ctrl-C' to exit");


%% Problem 1a: Prove SNC works with optimal sigma
    % Based on plots, sigma = 1e-8 balances both postfit and 3D RMS
fprintf("\n1a. Running LKF and Smoother with SNC for sigma = %.3e km/s^2\n", sigAccel)

    % Define optimal Q0
if measNoise
    Q0 = 1*diag(sigAccel^2*ones(1,3));
else
    Q0 = zeros(3,3);
end

    % Run LKF with SNC
LKFOpt = runLKF_SNC(X0, x0, P0, Q0, earthConst, stations, X_ref, t_ref, 1, false);

%% Problem 1b/d+e: Implement sequential filter smoothing algorithm and compare RMS values
fprintf("\n1b/d/e. Smoothing LKF SNC results\n")

    % Run Smoother algorithm
SmoothOpt = runSmoother(LKFOpt.LKFOut, X_ref, t_ref, true);

fprintf("\n\tLKF w/ SNC State RMS Errors:\n\t\tComponent-wise: [%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]\n\t\t3D: %.3e\n", LKFOpt.RMS_state_comp_LKF, LKFOpt.RMS_state_full_LKF);
fprintf("\n\tSmoother w/ SNC State RMS Errors:\n\t\tComponent-wise: [%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]\n\t\t3D: %.3e\n", SmoothOpt.RMS_state_comp_Smooth, SmoothOpt.RMS_state_full_Smooth);

%% Problem 1c: Smoothing without process noise
fprintf("\n1c. Running LKF w/ Smoother and Batch, both without Process noise\n")
    
    % LKF + Smoother
LKFNoProcess = runLKF_SNC(X0, x0, P0, 0*Q0, earthConst, stations, X_ref, t_ref, 1, true);
smoothNoProcess = runSmoother(LKFNoProcess.LKFOut, X_ref, t_ref, true);

    % Batch
BatchNoProcess = runBatch(X0, x0, P0, earthConst, stations, X_ref, t_ref, dt, opt, 1, []);





