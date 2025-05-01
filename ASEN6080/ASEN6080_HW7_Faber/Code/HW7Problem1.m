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
filename = "..\Data\HW7Data_J2J3.mat";

    % Only generate data once
if true
    generateTruthData_MuJ2J3(earthConst, orbital, x0Perturb, stations, filename, numOrbits, measNoise, dt);
end

load("..\Data\HW7Data_J2J3.mat"); % Source of X0, updated stations, X_ref, and t_ref

titleText = sprintf("Measurement data");
xLabel = "Time [sec]"; 
yLabel = ["\rho [km]", "\rhoDot [km/s]", "Elevation Angle [deg]"];
plotMeasurements(stations, titleText, xLabel, yLabel);

% input("Press 'Enter' to continue")

%% Problem 1a. Filter setup
    % Covariance matrices
Pxx0 = 1e4*eye(6);
Pcc0 = earthConst.J3^2;

    % Initial state deviation estimate
x0 = zeros(size(X0));

    % Initial sensitivity matrix
[~,~,~,~,vis] = processStations(stations);
X_s = stations(vis{1}).Xs(1,1:6); % Whichever station is first visible
[Hx, Hc] = MeasurementPartials_CP_J3_sc(X0, X_s);
S0 = -Pxx0*Hx'*stations(1).R{1}^-1*Hc;

    % Consider parameter setup
C = 0;
c = earthConst.J3 - C;

    % Plotting
plot = [true; true; false; false; false; false; true; true]; % Only plot residuals and state errors

%% Problem 1b/c. Run Consider Covariance filter
fprintf("\nRunning CFA on data\n")
% CFA_test = CFA(X0, stations, earthConst, Pxx0, Pcc0, x0, S0);

CFARun = runCFA(X0, x0, Pxx0, Pcc0, S0, earthConst, stations, X_ref, t_ref, plot);

%% Problem 1d/e. Map xHat_f and Pc,f to t0 and propagate again
fprintf("\nBackpropagating CFA solution\n\n")

    % Extract results from CFA run
xf = CFARun.CFAOut.xEst(:,end);
Pcf = CFARun.CFAOut.PcEst{end};
Pxcf = CFARun.CFAOut.PxcEst{end};
Psi = CFARun.CFAOut.Psi;

    % Accumulate final covariance
Pf = [
        Pcf, Pxcf;
        Pxcf', Pcc0
     ];

    % Propagate xf and Pcf back in time to t0
x0_new = CFARun.CFAOut.Phi_total{end}^-1*xf;
Pc0_new = (Psi{end}^-1)*Pf*(Psi{end}^-1)';

    % Propagate new initial state and covariance forward in time
X0_new = X0 + x0_new;
[t_new, X_new] = ode45(@(t,X)orbitEOM_MuJ2(t,X,earthConst.mu,earthConst.J2,earthConst.Ri), t_ref, X0_new, opt);

Pc_new = [];
for k = 1:length(Psi)
    P_part = Psi{k}*Pc0_new*Psi{k}';
    Pc_new = [Pc_new, {P_part}];
end

    % Calculate new state errors and find sigmas
if false
        % Create new state estimate vectors that align with the measured
        % covariance estimates
    X_new_CFA = [];
    X_ref_CFA = [];
    for k = 1:length(CFARun.t_CFA)
        X_new_CFA = [X_new_CFA; X_new(t_ref == CFARun.t_CFA(k),:)];
        X_ref_CFA = [X_ref_CFA; X_ref(t_ref == CFARun.t_CFA(k),:)];
    end
    t_CFA = CFARun.t_CFA;
else
    X_new_CFA = X_new;
    X_ref_CFA = X_ref;
    t_CFA = t_new;
end

        % Calculate errors
stateError_new = X_new_CFA - X_ref_CFA;

        % Find sigmas
sigma_new = [];
for k = 1:length(Pc_new)
    P = Pc_new{k};
        
    sigPart = [];
    for kk = 1:size(P,1)-1
        sigPart = [sigPart, sqrt(P(kk,kk))];
    end

    sigma_new = [sigma_new; sigPart];
end

    % Find new state error RMS
[RMS_state_comp_new, RMS_state_full_new] = calcStateErrorRMS(stateError_new);

fprintf("State Error Component-wise RMS, backpropagated: [%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]\n", RMS_state_comp_new);
fprintf("State Error 3D RMS, backpropagated: %.3e\n", RMS_state_full_new);

    % Plot new state errors with 2-sigma bounds
idx = 1;
boundLevel = 2; % Plot +/- boundLevel*sigma around state errors
titleText = "CFA Estimated State Error - backpropagated (X_{filt} - X_{ref})";
xLabel = "Time [sec]";
yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
          "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]"];

plotStateError(t_CFA(idx:end), stateError_new(idx:end,:), CFARun.t_CFA(idx:end), sigma_new(idx:end,:), boundLevel, titleText, xLabel, yLabel);


