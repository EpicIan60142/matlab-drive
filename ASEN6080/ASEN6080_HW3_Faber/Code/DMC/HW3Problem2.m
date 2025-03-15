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

    save("DMCData-new.mat", "LKFRuns", "EKFRuns", '-mat');
else
    load("DMCData.mat")
end

%% Problem 2b: Plot results vs. sigma
fprintf("\nPlotting DMC filter results vs. sigma\n")

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

%% Problem 2b: Choose optimal sigma
    % Based on plots, sigma = 1e-10 balances both postfit and 3D RMS
sigOpt = find(round(sigAccel,20) == 1e-9);

fprintf("\nPlotting state errors vs. time for sigma = %.3e km/s^2\n", sigAccel(sigOpt))

    % Set up optimal filter variables
sig_u = sigAccel(sigOpt);

Pw0 = diag([sig_u^2, sig_u^2, sig_u^2]);
P0_new = blkdiag(P0,Pw0);

    % Define process noise covariance matrix
Qu = diag([sig_u^2, sig_u^2, sig_u^2]);

    % Run filters
LKFOpt = runLKF_DMC(X0, x0, P0_new, B, Qu, earthConst, stations, X_ref, t_ref, 1, true);

numMeas = 100;
X0_EKF = LKFOpt.X_LKF(:,numMeas); P0_EKF = LKFOpt.LKFOut.PEst{numMeas}; t_start = LKFOpt.t_LKF(numMeas);
EKFOpt = runEKF_DMC(X0_EKF, P0_EKF, B, Qu, earthConst, stations, X_ref, t_ref, t_start, true);

%% Problem 2b: Plot J3 accels and estimated filter estimates vs. time
fprintf("\nPlotting true J3 accels vs. filter estimates\n")

    % Pull out estimated J3 accelerations
LKFOpt.J3Est = LKFOpt.X_LKF(7:9,:);
EKFOpt.J3Est = EKFOpt.X_EKF(7:9,:);

    % Calculate true J3 acceleration
mu = earthConst.mu;
J3 = earthConst.J3;
Ri = earthConst.Ri;
J3Func = @(x,y,z,r,Ri,mu,J3) ...
            [
                -(mu*x/(r^3))*((((Ri/r)^3)*J3*((35/2)*(z/r)^3-(15/2)*(z/r))));
                -(mu*y/(r^3))*((((Ri/r)^3)*J3*((35/2)*(z/r)^3-(15/2)*(z/r))));
                -(mu/(r^2))*((Ri/r)^3)*J3*((35/2)*(z/r)^4-15*(z/r)^2+(3/2));
            ];

J3Accel = [];
for k = 1:size(X_ref,1)
    x = X_ref(k,1); y = X_ref(k,2); z = X_ref(k,3);
    r = sqrt(x^2 + y^2 + z^2);

    J3Accel = [J3Accel, J3Func(x,y,z,r,Ri,mu,J3)];
end

    % Find difference of true and estimated J3
LKFOpt.J3Diff = [];
EKFOpt.J3Diff = [];
for k = 1:length(t_ref)
        % LKF
    idx = LKFOpt.t_LKF == t_ref(k);
    if ~isempty(idx)
        LKFOpt.J3Diff = [LKFOpt.J3Diff, J3Accel(:,k) - LKFOpt.J3Est(:,idx)];
    end
        
        % EKF
    idx = EKFOpt.t_EKF == t_ref(k);
    if ~isempty(idx)
        EKFOpt.J3Diff = [EKFOpt.J3Diff, J3Accel(:,k) - EKFOpt.J3Est(:,idx)];
    end
end

    % Plot estimated and true J3 for each filter
figure; tl = tiledlayout(3,1); ax = [];
title(tl,"J3 Accelerations vs. time")
titles = ["J_{3,X} vs. Time", "J_{3,Y} vs. Time", "J_{3,Z} vs. Time"]; 
yLabels = ["J_{3,X}", "J_{3,Y}", "J_{3,Z}"];
for k = 1:3
    nt = nexttile; ax = [ax; nt];
        title(titles(k))
        hold on; grid on;
        plot(LKFOpt.t_LKF, LKFOpt.J3Est(k,:), 'b.');
        plot(EKFOpt.t_EKF, EKFOpt.J3Est(k,:), 'r.');
        plot(t_ref, J3Accel(k,:), 'k-');
        xlabel("Time [sec]"); ylabel(yLabels(k));
end
legend("LKF Estimated", "EKF Estimated", "True J_3", 'Location', 'northeastoutside')
linkaxes(ax,'x')

    % Plot J3 differences with covariance bounds
        % Find sigma bounds for LKF estimated unmodeled accelerations
sigma_LKF = [];
PEst_LKF = LKFOpt.LKFOut.PEst;
for k = 1:length(PEst_LKF)
    P = PEst_LKF{k};
    P = P(7:9, 7:9);

    sigPart = [];
    for kk = 1:size(P,1)
        sigPart = [sigPart, sqrt(P(kk,kk))];
    end

    sigma_LKF = [sigma_LKF; sigPart];
end
   
        % find sigma bounds for EKF estimated unmodeled accelerations
sigma_EKF = [];
PEst_EKF = EKFOpt.EKFOut.PEst;
for k = 1:length(PEst_EKF)
    P = PEst_EKF{k};
    P = P(7:9, 7:9);

    sigPart = [];
    for kk = 1:size(P,1)
        sigPart = [sigPart, sqrt(P(kk,kk))];
    end

    sigma_EKF = [sigma_EKF; sigPart];
end

    % Plot errors with bounds
titleText = "LKF Estimated Unmodeled J_3 Error (J_{3,true} - J_{3,LKF})";
xLabel = "Time [sec]";
yLabel = ["J_{3,X} error [km/s^2]", "J_{3,Y} error [km/s^2]", "J_{3,Z} error [km/s^2]"];
plotStateError(LKFOpt.t_LKF, LKFOpt.J3Diff', LKFOpt.t_LKF, sigma_LKF, 3, titleText, xLabel, yLabel);

titleText = "EKF Estimated Unmodeled J_3 Error (J_{3,true} - J_{3,EKF})";
xLabel = "Time [sec]";
yLabel = ["J_{3,X} error [km/s^2]", "J_{3,Y} error [km/s^2]", "J_{3,Z} error [km/s^2]"];
plotStateError(EKFOpt.t_EKF, EKFOpt.J3Diff', EKFOpt.t_EKF, sigma_EKF, 3, titleText, xLabel, yLabel);

return

%% Problem 2b: Different value of tau
tau_x = T/30; tau_y = T/30; tau_z = T/600;
B = diag([tau_x^-1, tau_y^-1, tau_z^-1]);

fprintf("Running filters with tau_x = %.3e, tau_y = %.3e, tau_z = %.3e", tau_x, tau_y, tau_z);

    % Set up optimal filter variables
sig_u = sigAccel(sigOpt);

Pw0 = diag([sig_u^2, sig_u^2, sig_u^2]);
P0_new = blkdiag(P0,Pw0);

    % Define process noise covariance matrix
Qu = diag([sig_u^2, sig_u^2, sig_u^2]);

    % Run filters
LKFNewTau = runLKF_DMC(X0, x0, P0_new, B, Qu, earthConst, stations, X_ref, t_ref, 5, true);

numMeas = 157;
X0_EKF = LKFNewTau.X_LKF(:,numMeas); P0_EKF = LKFNewTau.LKFOut.PEst{numMeas}; t_start = LKFNewTau.t_LKF(numMeas);
EKFNewTau = runEKF_DMC(X0_EKF, P0_EKF, B, Qu, earthConst, stations, X_ref, t_ref, t_start, true);

    % Pull out estimated J3 accelerations
LKFNewTau.J3Est = LKFNewTau.X_LKF(7:9,:);
EKFNewTau.J3Est = EKFNewTau.X_EKF(7:9,:);

    % Find difference of true and estimated J3
LKFNewTau.J3Diff = [];
EKFNewTau.J3Diff = [];
for k = 1:length(t_ref)
        % LKF
    idx = LKFNewTau.t_LKF == t_ref(k);
    if ~isempty(idx)
        LKFNewTau.J3Diff = [LKFNewTau.J3Diff, J3Accel(:,k) - LKFNewTau.J3Est(:,idx)];
    end
        
        % EKF
    idx = EKFNewTau.t_EKF == t_ref(k);
    if ~isempty(idx)
        EKFNewTau.J3Diff = [EKFNewTau.J3Diff, J3Accel(:,k) - EKFNewTau.J3Est(:,idx)];
    end
end

    % Plot estimated and true J3 for each filter
figure; tl = tiledlayout(3,1); ax = [];
title(tl,"J3 Accelerations vs. time")
titles = ["J_{3,X} vs. Time", "J_{3,Y} vs. Time", "J_{3,Z} vs. Time"]; 
yLabels = ["J_{3,X}", "J_{3,Y}", "J_{3,Z}"];
for k = 1:3
    nt = nexttile; ax = [ax; nt];
        title(titles(k))
        hold on; grid on;
        plot(LKFNewTau.t_LKF, LKFNewTau.J3Est(k,:), 'b.');
        plot(EKFNewTau.t_EKF, EKFNewTau.J3Est(k,:), 'r.');
        plot(t_ref, J3Accel(k,:), 'k-');
        xlabel("Time [sec]"); ylabel(yLabels(k));
end
legend("LKF Estimated", "EKF Estimated", "True J_3", 'Location', 'bestoutside')
linkaxes(ax,'x')

    % Plot errors with bounds
titleText = "LKF Estimated Unmodeled J_3 Error (J_{3,true} - J_{3,LKF})";
xLabel = "Time [sec]";
yLabel = ["J_{3,X} error [km/s^2]", "J_{3,Y} error [km/s^2]", "J_{3,Z} error [km/s^2]"];
plotStateError(LKFNewTau.t_LKF, LKFNewTau.J3Diff', LKFNewTau.t_LKF, sigma_LKF, 3, titleText, xLabel, yLabel);

titleText = "EKF Estimated Unmodeled J_3 Error (J_{3,true} - J_{3,EKF})";
xLabel = "Time [sec]";
yLabel = ["J_{3,X} error [km/s^2]", "J_{3,Y} error [km/s^2]", "J_{3,Z} error [km/s^2]"];
plotStateError(EKFNewTau.t_EKF, EKFNewTau.J3Diff', EKFNewTau.t_EKF, sigma_EKF, 3, titleText, xLabel, yLabel);


