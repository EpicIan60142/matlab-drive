%% ASEN 6080 HW 3 Problem 2 Main Script
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

%% Problem 2b: Filter setup
sigR = 1; % km
sigV = 1e-3; % km/s

eta = zeros(3,1);
X0 = [X0; eta]; % Append a priori w states
x0 = 0*0.5*[1 1 1 1e-3 1e-3 1e-3 0 0 0]';
P0 = diag([sigR^2, sigR^2, sigR^2, sigV^2, sigV^2, sigV^2]);

T = 9520; % Approximate orbit period

tau_x = T/30; tau_y = T/30; tau_z = T/30;
B = diag([tau_x^-1, tau_y^-1, tau_z^-1]);

sigAccel = logspace(-15,-2,14); % Create test sigmas from 1e-15 to 1e-2 km/s^2

numMeas = 50; % Number of LKF measurements to initialize EKF with

input("Press Enter to continue, 'Ctrl-C' to exit");

%% Problem 2b: Test values of sigma
if false
    LKFRuns = [];
    EKFRuns = [];
    for k = 1:length(sigAccel)
            % Pull out new sigma to test and update initial covariance
        sig_u = sigAccel(k);

        Pw0 = diag([sig_u^2, sig_u^2, sig_u^2]);
        P0_new = blkdiag(P0,Pw0);

        fprintf("\nRunning with sigmaAccel = %.2e:\n", sig_u);

            % Define process noise covariance matrix
        Qu = diag([sig_u^2, sig_u^2, sig_u^2]);
    
            % Run filters
        LKFRuns = [LKFRuns; runLKF_DMC(X0, x0, P0_new, B, Qu, earthConst, stations, X_ref, t_ref, 5, false)];
    
        LKFRun = LKFRuns(k);
        X0_EKF = LKFRun.X_LKF(:,numMeas); P0_EKF = LKFRun.LKFOut.PEst{numMeas}; t_start = LKFRun.t_LKF(numMeas);
        EKFRuns = [EKFRuns; runEKF_DMC(X0_EKF, P0_EKF, B, Qu, earthConst, stations, X_ref, t_ref, t_start, false)];
    end

    save("DMCData.mat", "LKFRuns", "EKFRuns", '-mat');
else
    load("DMCData.mat")
end

%% Problem 2b: Plot results vs. sigma
fprintf("\nPlotting SNC filter results vs. sigma\n")

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
xlabel("\sigma_{accel} [km/s^2]"); ylabel("Normalized Postfit RMS")
legend("LKF", "EKF", 'location', 'eastoutside');

figure; 
loglog(sigAccel, LKF_3DRMS, 'b-');
hold on; grid on;
loglog(sigAccel, EKF_3DRMS, 'r-');
title("3-D RMS vs. \sigma_{accel}")
xlabel("\sigma_{accel} [km/s^2]"); ylabel("3-D RMS")
legend("LKF", "EKF", 'location', 'eastoutside');
drawnow;







