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

%% Load Truth Data
load("..\Data\HW3Data.mat");

titleText = sprintf("Measurement data");
xLabel = "Time [sec]"; 
yLabel = ["\rho [km]", "\rhoDot [km/s]", "Elevation Angle [deg]"];
plotMeasurements(stations, titleText, xLabel, yLabel);

%% Problem 1b: Filter setup
sigR = 1; % km
sigV = 1e-3; % km/s

x0 = 0*0.5*[1 1 1 1e-3 1e-3 1e-3]';
P0 = diag([sigR^2, sigR^2, sigR^2, sigV^2, sigV^2, sigV^2]);

sigAccel = logspace(-15,-2,14); % Create test sigmas from 1e-15 to 1e-2 m/s^2

numMeas = 50; % Number of LKF measurements to initialize the EKF with

input("Press Enter to continue, 'Ctrl-C' to exit");

%% Problem 1b: Run filters once per sigma 
if false
    LKFRuns = [];
    EKFRuns = [];
    for k = 1:length(sigAccel)
            % Pull out next sigma
        sigma = sigAccel(k);
    
        fprintf("\nRunning with sigmaAccel = %.2e:\n", sigma);
    
            % Define Q0
        Q0 = diag(sigma^2*ones(1,3));
    
            % Run LKF
        LKFRuns = [LKFRuns; runLKF_SNC(X0, x0, P0, Q0, earthConst, stations, X_ref, t_ref, 5, false)];
    
            % Run EKF
        LKFRun = LKFRuns(k);
        X0_EKF = LKFRun.X_LKF(:,numMeas); P0_EKF = LKFRun.LKFOut.PEst{numMeas}; t_start = LKFRun.t_LKF(numMeas);
        EKFRuns = [EKFRuns; runEKF_SNC(X0_EKF,P0_EKF,Q0,earthConst, stations, X_ref, t_ref, t_start, false)];
    end
    
    save("SNCData.mat", "LKFRuns", "EKFRuns", '-mat');
else
    load("SNCData.mat");
end
%% Problem 1b: Plot results vs. sigma
    % Compile RMS
LKF_PostFitRMS = []; LKF_3DRMS = [];
EKF_PostFitRMS = []; EKF_3DRMS = [];
for k = 1:length(LKFRuns)
    LKF_PostFitRMS = [LKF_PostFitRMS; LKFRuns(k).RMS_postfit_LKF(end)];
    LKF_3DRMS = [LKF_3DRMS; LKFRuns(k).RMS_state_full_LKF(end)];
    EKF_PostFitRMS = [EKF_PostFitRMS; EKFRuns(k).RMS_postfit_EKF(end)];
    EKF_3DRMS = [EKF_3DRMS; EKFRuns(k).RMS_state_full_EKF(end)];
end

    % Plot RMS vs. sigma
figure; 
loglog(sigAccel, LKF_PostFitRMS, 'b-');
hold on; grid on;
loglog(sigAccel, EKF_PostFitRMS, 'r-');
title("Normalized Post-Fit RMS vs. \sigma_{accel}")
xlabel("\sigma_{accel} [m/s^2]"); ylabel("Normalized Postfit RMS")
legend("LKF", "EKF", 'location', 'eastoutside');

figure; 
loglog(sigAccel, LKF_3DRMS, 'b-');
hold on; grid on;
loglog(sigAccel, EKF_3DRMS, 'r-');
title("3-D RMS vs. \sigma_{accel}")
xlabel("\sigma_{accel} [m/s^2]"); ylabel("3-D RMS")
legend("LKF", "EKF", 'location', 'eastoutside');
drawnow;

%% Problem 1b: Choose optimal sigma and show state errors
    % Based on plots, sigma = 1e-8 balances both postfit and 3D RMS
sigOpt = find(sigAccel == 1e-8); 

    % Define optimal Q0
Q0 = diag(sigAccel(sigOpt)^2*ones(1,3));

    % Run filters
LKFOpt = runLKF_SNC(X0, x0, P0, Q0, earthConst, stations, X_ref, t_ref, 5, true);

numMeas = 50; % Initialize with 50 LKF measurements
X0_EKF = LKFOpt.X_LKF(:,numMeas); P0_EKF = LKFOpt.LKFOut.PEst{numMeas}; t_start = LKFOpt.t_LKF(numMeas);
EKFOpt = runEKF_SNC(X0_EKF, P0_EKF, Q0, earthConst, stations, X_ref, t_ref, t_start, true);

