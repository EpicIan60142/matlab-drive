%% ASEN 6080 Project 2 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
    % Path logistics
addpath("..\")
addpath(genpath("..\..\Utilities\"));

    % Celestial body constants
pConst = getPlanetConst();

    % Spacecraft constants
scConst = getSCConst();

    % Initialize stations struct
stations = makeStations(pConst, 2);

    % ode45 options
opt = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

%% Part 1: Dynamics Implementation
    % Set up initial state and get initial accelerations
X0 = [scConst.X0_cart; scConst.C_R];
aTest = orbitEOM_MuSunSRP(0, X0, pConst, scConst);

    % Load truth trajectory to verify dynamics are working
truthTraj = load("..\Data\Project2_Prob2_truth_traj_50days.mat");

    % Integrate X0 to verify EOM function
tspan = [truthTraj.Tt_50; (truthTraj.Tt_50(end)+1000:1000:279*24*60*60)']; % Add time up to 200 days in sec
[t_test, X_test] = ode45(@(t,X)orbitEOM_MuSunSRP(t,X,pConst,scConst), tspan, truthTraj.Xt_50(1,1:7)', opt);

    % Plot error between modeled and true trajectory
idx = 1:length(truthTraj.Tt_50); % Only plot data when truth data exists

figure; tl = tiledlayout(2,3);
title(tl, "Part 1: Error between modeled and true trajectory")
nexttile;
    hold on; grid on;
    plot(t_test(idx), X_test(idx,1) - truthTraj.Xt_50(idx,1), 'b-');
    xlabel("Time [sec]"); ylabel("\DeltaX [km]")
nexttile;
    hold on; grid on;
    plot(t_test(idx), X_test(idx,2) - truthTraj.Xt_50(idx,2), 'b-');
    xlabel("Time [sec]"); ylabel("\DeltaY [km]")
nexttile;
    hold on; grid on;
    plot(t_test(idx), X_test(idx,3) - truthTraj.Xt_50(idx,3), 'b-');
    xlabel("Time [sec]"); ylabel("\DeltaZ [km]")
nexttile;
    hold on; grid on;
    plot(t_test(idx), X_test(idx,4) - truthTraj.Xt_50(idx,4), 'b-');
    xlabel("Time [sec]"); ylabel("\DeltaXdot [km/s]")
nexttile;
    hold on; grid on;
    plot(t_test(idx), X_test(idx,5) - truthTraj.Xt_50(idx,5), 'b-');
    xlabel("Time [sec]"); ylabel("\DeltaYdot [km/s]")
nexttile;
    hold on; grid on;
    plot(t_test(idx), X_test(idx,6) - truthTraj.Xt_50(idx,6), 'b-');
    xlabel("Time [sec]"); ylabel("\DeltaZdot [km/s]")

    % Plot initial situation
[R_Earth, ~, ~, ~] = Ephem(pConst.initEpoch, 3, 'EME2000');
R_E_t = []; R_E_t_2 = [];
for k = 1:length(t_test)%truthTraj.Tt_50)
    % [R_E_t(:,k), ~, ~, ~] = Ephem(pConst.initEpoch + truthTraj.Tt_50(k)/(24*60*60), 3, 'EME2000'); % Earth orbit over this scenario
    [R_E_t(:,k), ~, ~, ~] = Ephem(pConst.initEpoch + t_test(k)/(24*60*60), 3, 'EME2000'); % Earth orbit over this scenario
end

for k = 1:365
    [R_E_t_2(:,k), ~, ~, ~] = Ephem(pConst.initEpoch + k, 3, 'EME2000'); % Earth full orbit
end

figure;
title("Part 1: Truth orbit comparison until 200 days - Sun origin frame")
hold on; grid on; axis equal;
plot3(0,0,0,'y.','MarkerSize', 50);
plot3(X_test(1,1)+R_Earth(1), X_test(1,2)+R_Earth(2), X_test(1,3)+R_Earth(3),'m.','MarkerSize',30);
plot3(truthTraj.Xt_50(:,1)+R_E_t(1,1:length(truthTraj.Tt_50))', truthTraj.Xt_50(:,2)+R_E_t(2,1:length(truthTraj.Tt_50))', truthTraj.Xt_50(:,3)+R_E_t(3,1:length(truthTraj.Tt_50))', 'm-.', 'LineWidth', 3)
plot3(X_test(:,1)+R_E_t(1,:)', X_test(:,2)+R_E_t(2,:)', X_test(:,3)+R_E_t(3,:)', 'm:', 'LineWidth', 2)
plot3(R_Earth(1), R_Earth(2), R_Earth(3),'b.','MarkerSize',40);
plot3(R_E_t(1,:), R_E_t(2,:), R_E_t(3,:), 'b-.', 'LineWidth', 2)
plot3(R_E_t_2(1,:), R_E_t_2(2,:), R_E_t_2(3,:), 'k:')
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]")
legend("Sun", "SC", "SC Truth Orbit", "SC Modeled Orbit", "Earth", "Earth Orbit")
view([30 35])

figure;
title("Part 1: Truth orbit comparison until 200 days - Earth origin frame")
hold on; grid on; axis equal;
plot3(0,0,0,'b.','MarkerSize', 50);
plot3(X_test(1,1), X_test(1,2), X_test(1,3),'m.','MarkerSize',30);
plot3(truthTraj.Xt_50(:,1), truthTraj.Xt_50(:,2), truthTraj.Xt_50(:,3), 'm-.', 'LineWidth', 3)
plot3(X_test(:,1), X_test(:,2), X_test(:,3), 'm:', 'LineWidth', 2)
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]")
legend("Earth", "SC", "SC Truth Orbit", "SC Modeled Orbit")
view([30 35])

%% Part 1: DSN Measurement Processing
data2a = readmatrix("..\Data\Project2a_Obs.txt");
[tMeas, stations] = readObsData(stations, data2a);

    % Propagate station states and generate clean measurements
for k = 1:length(tMeas)
    dTheta = pConst.wEarth*t_test(k);
    for kk = 1:length(stations)
        r = rotZ(dTheta)*stations(kk).X0;
        v = cross([0;0;pConst.wEarth], r);

        stations(kk).Xs = [stations(kk).Xs; [r', v', tMeas(k)]];        
    end
end

stations_nom = makeStations(pConst, 2);
    % Propagate station states and generate clean measurements
for k = 1:length(t_test)
    dTheta = pConst.wEarth*t_test(k);
    for kk = 1:length(stations_nom)
        r = rotZ(dTheta)*stations_nom(kk).X0;
        v = cross([0;0;pConst.wEarth], r);

        stations_nom(kk).Xs = [stations_nom(kk).Xs; [r', v', t_test(k)]];
        
        y = generateRngRngRate(X_test(k,:), stations_nom(kk).Xs(k,:), stations_nom(kk).elMask);

        if ~isnan(y)
            stations_nom(kk).rho = [stations_nom(kk).rho; y(1)];
            stations_nom(kk).rhoDot = [stations_nom(kk).rhoDot; y(2)];
            stations_nom(kk).elAngle = [stations_nom(kk).elAngle; y(3)];
            stations_nom(kk).tMeas = [stations_nom(kk).tMeas; t_test(k)];
            if true
                R = diag([stations_nom(kk).sigRho^2, stations_nom(kk).sigRhoDot^2]);
                stations_nom(kk).R = [stations_nom(kk).R; {R}];
            else
                R = zeros(2,2);
                stations_nom(kk).R = [stations_nom(kk).R; {R}];
            end
        end
    end
end

titleText = sprintf("Part 1: Provided Measurement Data");
xLabel = "Time [sec]"; 
yLabel = ["\rho [km]", "\rhoDot [km/s]"];
plotMeasurements(stations, titleText, xLabel, yLabel);

titleText = sprintf("Part 1: Propagated Measurement Data");
xLabel = "Time [sec]"; 
yLabel = ["\rho [km]", "\rhoDot [km/s]"];%, "Elevation Angle [deg]"];
plotMeasurements(stations_nom, titleText, xLabel, yLabel);

%% Part 1: Bplane implementation
    % Set up ODE events for sphere of influence
opt.Events = @(t,X)SOICheck(t,X,pConst);

    % Define function handle for integration
DynFunc = @(t,XPhi)STMEOM_MuSunSRP(t,XPhi,pConst,scConst);

    % Integrate to 3 SOI
XPhi_0 = [truthTraj.Xt_50(1,1:7)'; reshape(eye(7), 49, 1)];
P0 = 1e-11*eye(7); % Test covariance - not physical
tspan_3SOI = [truthTraj.Tt_50; (truthTraj.Tt_50(end)+1000:1000:300*24*60*60)']; % Add time up to 300 days in sec
[t_3SOI, XPhi_3SOI] = ode45(@(t,XPhi)DynFunc(t,XPhi), tspan_3SOI, XPhi_0, opt);
Phi_3SOI = reshape(XPhi_3SOI(end, 8:end),7,7);
P_3SOI = Phi_3SOI*P0*Phi_3SOI';

    % Calculate Bplane without SOI checker
[BdotR, BdotT, X_crossing, P_BPlane, STR2ECI, XPhi_BPlane, t_BPlane] = calcBPlane(XPhi_3SOI(end,:)', t_3SOI(end), P_3SOI, pConst, DynFunc, odeset('RelTol',1e-13,'AbsTol',1e-13));

    % Plot BdotR and BdotT location + uncertainty ellipse
titleText = sprintf("B Plane Crossing Estimate");
xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]";
plotBPlane(BdotR, BdotT, X_crossing, P_BPlane, STR2ECI, pConst, 3, titleText, xLabel, yLabel, zLabel, 420);

figure(420)
hold on;
idx = length(t_BPlane); offset = 1000;
orbit = plot3(XPhi_BPlane(idx-offset:end,1), XPhi_BPlane(idx-offset:end,2), XPhi_BPlane(idx-offset:end,3), 'm--', 'DisplayName', "Modeled Orbit");

return

%% Part 2: Estimate State with Known Target and Models
    % Set initial state estimate and deviation
% X0 = [scConst.X0_cart; scConst.C_R];
X0 = truthTraj.Xt_50(end,1:7)';
x0 = zeros(size(X0));

    % Set initial covariances
sigR = 100; % km
sigV = 0.1; % km/s
sigCR = 0.1;
P0 = diag([sigR^2, sigR^2, sigR^2, sigV^2, sigV^2, sigV^2, sigCR^2]);

    % Reset ODE45 options - won't need SOI checker
opt = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

    % Remake and populate stations structure
        % Make structure
stations = makeStations(pConst, 2);
        % Read and populate observations
data2a = readmatrix('..\Data\Project2a_Obs.txt');
[tMeas, stations] = readObsData(stations, data2a);

        % Propagate station states
for k = 1:length(tMeas)
    dTheta = pConst.wEarth*t_test(k);
    for kk = 1:length(stations)
        r = rotZ(dTheta)*stations(kk).X0;
        v = cross([0;0;pConst.wEarth], r);

        stations(kk).Xs = [stations(kk).Xs; [r', v', tMeas(k)]];        
    end
end

    % Create new truth data for first 50 days
kEnd = find(tMeas <= 50*24*60*60, 1, 'last'); 
[t_50, X_50] = ode45(@(t,X)orbitEOM_MuSunSRP(t,X,pConst,scConst), tMeas(1:kEnd), truthTraj.Xt_50(1,1:7),opt);

    % Run UKF on 50 days of data
tSpan = [0, 50*24*60*60]; % 0 to 50 days in seconds
alpha = 1; beta = 2;
plot = [true; true; false; false; false; false; true; true]; % Only plot residuals and state errors
UKFRun = runUKF(X0, P0, zeros(3,3), tSpan, pConst, scConst, stations, X_50, t_50, alpha, beta, opt, plot);


