%% ASEN 6080 HW 2 Problem 1 Main Script
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

%% Make Truth Data
x0Perturb = 1*0.5*[1 1 1 1e-3 1e-3 1e-3]';
numOrbits = 15;
dt = 10; % sec
filename = "HW2Problem1Data-new.mat";

    % Only generate data once
if true
    generateTruthData_MuJ2(earthConst, orbital, 0*x0Perturb, stations, filename, numOrbits, dt);
end

load("HW2Problem1Data-new.mat");

titleText = sprintf("Measurement data");
xLabel = "Time [sec]"; 
yLabel = ["\rho [km]", "\rhoDot [km/s]", "Elevation Angle [deg]"];
plotMeasurements(stations, titleText, xLabel, yLabel);

%% Problem 1a: Filter setup
sigR = 1; % km
sigV = 1e-3; % km/s

x0 = 1*0.5*[1 1 1 1e-3 1e-3 1e-3]';
P0 = diag([sigR^2, sigR^2, sigR^2, sigV^2, sigV^2, sigV^2]);

input("Press Enter to continue, 'Ctrl-C' to exit");

%% Problem 1b: Batch Filter
fprintf("1b. Batch Filter\n\n")

batchRun = runBatch(X0+x0, x0, P0, earthConst, stations, X_ref, tspan, dt, opt);

%% Problem 1b: CKF/LKF
fprintf("\n1b. LKF\n\n")

LKFRun = runLKF(X0+x0, x0, P0, earthConst, stations, X_ref, t_ref);

X_LKF = LKFRun.X_LKF;
PEst_LKF = LKFRun.LKFOut.PEst;
t_LKF = LKFRun.t_LKF;

%% Problem 1b: EKF
fprintf("\n1b. EKF\n\n")

    % Initialize EKF using some number of LKF measurements
numMeas = 157;
X0_EKF = X_LKF(:,numMeas);
x0_EKF = LKFRun.LKFOut.xEst(:,numMeas);
P0_EKF = PEst_LKF{numMeas};
t_EKF_start = t_LKF(numMeas);

EKFRun = runEKF(X0_EKF+x0_EKF, P0_EKF, earthConst, stations, X_ref, t_ref, t_EKF_start);

%% Problem 1d: Increase x0
x0_d = 100000*x0;

fprintf("\n1d. Increased x0\n\n")

batchRun_d = runBatch(X0+x0, x0_d, P0, earthConst, stations, X_ref, tspan, dt, opt);

LKFRun_d = runLKF(X0+x0, x0_d, P0, earthConst, stations, X_ref, t_ref);

numMeas = 157;
X0_EKF = LKFRun_d.X_LKF(:,numMeas);
x0_EKF = LKFRun_d.LKFOut.xEst(:,numMeas);
P0_EKF = LKFRun_d.LKFOut.PEst{numMeas};
t_EKF_start = LKFRun_d.t_LKF(numMeas);
EKFRun_d = runEKF(X0_EKF, P0_EKF, earthConst, stations, X_ref, t_ref, t_EKF_start);

%% Problem 1e: Increase P0
sigR_e = 1000; % km
sigV_e = 1; % km/s

P0_e = diag([sigR_e^2, sigR_e^2, sigR_e^2, sigV_e^2, sigV_e^2, sigV_e^2]);

fprintf("\n1e. Increased P0\n\n")

batchRun_e = runBatch(X0+x0, x0_d, P0_e, earthConst, stations, X_ref, tspan, dt, opt);

LKFRun_e = runLKF(X0+x0, x0_d, P0_e, earthConst, stations, X_ref, t_ref);

numMeas = 157;
X0_EKF = LKFRun_e.X_LKF(:,numMeas);
x0_EKF = LKFRun_e.LKFOut.xEst(:,numMeas);
P0_EKF = LKFRun_e.LKFOut.PEst{numMeas};
t_EKF_start = LKFRun_e.t_LKF(numMeas);
EKFRun_e = runEKF(X0_EKF, P0_EKF, earthConst, stations, X_ref, t_ref, t_EKF_start);

%% Problem 1f: Process half the measurements
stations_f = stations;
for k = 1:length(stations)
    stations_f(k).tMeas = stations(k).tMeas(1:end/2);
    stations_f(k).rho = stations(k).rho(1:end/2);
    stations_f(k).rhoDot = stations(k).rhoDot(1:end/2);
    stations_f(k).elAngle = stations(k).elAngle(1:end/2);
    stations_f(k).R = stations(k).R(1:end/2);
end

fprintf("\n1f. Half measurements\n\n")

batchRun_f = runBatch(X0+x0, x0, P0, earthConst, stations_f, X_ref, tspan, dt, opt);

LKFRun_f = runLKF(X0+x0, x0, P0, earthConst, stations_f, X_ref, t_ref);

numMeas = 157;
X0_EKF = LKFRun_f.X_LKF(:,numMeas);
x0_EKF = LKFRun_f.LKFOut.xEst(:,numMeas);
P0_EKF = LKFRun_f.LKFOut.PEst{numMeas};
t_EKF_start = LKFRun_f.t_LKF(numMeas);
EKFRun_f = runEKF(X0_EKF+x0_EKF, P0_EKF, earthConst, stations_f, X_ref, t_ref, t_EKF_start);

