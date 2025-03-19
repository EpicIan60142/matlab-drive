%% ASEN 6080 HW 6 Problem 1 Main Script
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
sigR = 1; % km
sigV = 1e-3; % km/s

P0 = diag([sigR^2, sigR^2, sigR^2, sigV^2, sigV^2, sigV^2]);

sigAccel_UKF = 1e-7; % Test value based on intuition, seems to work well!
sigAccel_EKF = 1e-8; % Optimal SNC sigma from HW 3

Q0_UKF_init = diag(sigAccel_UKF^2*ones(1,3));
Q0_EKF_init = diag(sigAccel_EKF^2*ones(1,3));

plot = [true; true; false; false; false; false; true; true]; % Only plot residuals and state errors

%% Problem 1b/c. UKF test cases
    % alpha = 1, beta = 2, no process noise
Q0_UKF = 0*Q0_UKF_init;
alpha_UKF = 1;
beta_UKF = 2;

fprintf("\nProblem 1b. UKF with alpha = %.4f, beta = %.0f, no process noise\n", alpha_UKF, beta_UKF);

UKFCase1 = runUKF(X0,P0,Q0_UKF,earthConst,stations,X_ref,t_ref,alpha_UKF,beta_UKF,false,plot);

    % Run EKF for comparison
Q0_EKF = 0*Q0_EKF_init;
LKFRun = runLKF_SNC(X0,zeros(size(X0)),P0,Q0_EKF,earthConst,stations,X_ref,t_ref,1,false);
numMeas = 50;
X0_EKF = LKFRun.X_LKF(:,numMeas); P0_EKF = LKFRun.LKFOut.PEst{numMeas}; t_start = LKFRun.t_LKF(numMeas);
EKFNoProcNoise = runEKF_SNC(X0_EKF,P0_EKF,Q0_EKF,earthConst,stations,X_ref,t_ref,t_start,plot);

input("Press 'Enter' to continue")

    % Add process noise
Q0_UKF = Q0_UKF_init;

fprintf("\nProblem 1b. UKF with alpha = %.4f, beta = %.0f, with process noise\n", alpha_UKF, beta_UKF);

UKFCase2 = runUKF(X0,P0,Q0_UKF,earthConst,stations,X_ref,t_ref,alpha_UKF,beta_UKF,false,plot);

    % Run EKF for comparison
Q0_EKF = Q0_EKF_init;
LKFRun = runLKF_SNC(X0,zeros(size(X0)),P0,Q0_EKF,earthConst,stations,X_ref,t_ref,1,false);
numMeas = 50;
X0_EKF = LKFRun.X_LKF(:,numMeas); P0_EKF = LKFRun.LKFOut.PEst{numMeas}; t_start = LKFRun.t_LKF(numMeas);
EKFProcNoise = runEKF_SNC(X0_EKF,P0_EKF,Q0_EKF,earthConst,stations,X_ref,t_ref,t_start,plot);

input("Press 'Enter' to continue")

    % Set alpha = 1e-4, all else the same
alpha_UKF = 1e-4;

fprintf("\nProblem 1b. UKF with alpha = %.4f, beta = %.0f, with process noise\n", alpha_UKF, beta_UKF);

UKFCase3 = runUKF(X0,P0,Q0_UKF,earthConst,stations,X_ref,t_ref,alpha_UKF,beta_UKF,false,plot);

input("Press 'Enter' to continue")

%% Problem 1d. Investigate how robust UKF is to large state errors compared to EKF
fprintf("\nProblem 1d. Investigate performance of UKF vs. EKF for large initial state errors\n")

X0_UKF = 2*X0;

alpha_UKF = 1;

UKFLargeError = runUKF(X0_UKF,P0,Q0_UKF,earthConst,stations,X_ref,t_ref,alpha_UKF,beta_UKF,false,plot);

    % Run EKF for comparison
EKFLargeError = runEKF_SNC(X0_UKF,P0,Q0_EKF,earthConst,stations,X_ref,t_ref,0,plot);

input("Press 'Enter' to continue")

%% Problem 1e. Include J3 in filter dynamics
fprintf("\nProblem 1e. Add J3 to filter dynamics\n")
X0_UKF = X0;

UKFJ3 = runUKF(X0_UKF,P0,Q0_UKF,earthConst,stations,X_ref,t_ref,alpha_UKF,beta_UKF,true,plot);




