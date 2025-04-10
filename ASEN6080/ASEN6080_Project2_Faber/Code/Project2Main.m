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
% X0 = truthTraj.Xt_50(1,1:7)';

    % Integrate X0 to verify EOM function
tspan = [truthTraj.Tt_50; (truthTraj.Tt_50(end)+1000:1000:200*24*60*60)']; % Add time up to 200 days in sec
[t_test, X_test] = ode45(@(t,X)orbitEOM_MuSunSRP(t,X,pConst,scConst), tspan, truthTraj.Xt_50(1,1:7)', opt);

    % Plot error between modeled and true trajectory
idx = 1:length(truthTraj.Tt_50); % Only plot data when truth data exists

figure; tl = tiledlayout(2,3);
title(tl, "Error between modeled and true trajectory")
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
for k = 1:length(truthTraj.Tt_50)
    [R_E_t(:,k), ~, ~, ~] = Ephem(pConst.initEpoch + truthTraj.Tt_50(k)/(24*60*60), 3, 'EME2000'); % Earth orbit over this scenario
end

for k = 1:365
    [R_E_t_2(:,k), ~, ~, ~] = Ephem(pConst.initEpoch + k, 3, 'EME2000'); % Earth full orbit
end

figure;
title("Initial time situation")
hold on; grid on; axis equal;
plot3(0,0,0,'y.','MarkerSize', 50);
plot3(X_test(1,1)+R_Earth(1), X_test(1,2)+R_Earth(2), X_test(1,3)+R_Earth(3),'m.','MarkerSize',30);
plot3(truthTraj.Xt_50(:,1)+R_Earth(1), truthTraj.Xt_50(:,2)+R_Earth(2), truthTraj.Xt_50(:,3)+R_Earth(3), 'm-.', 'LineWidth', 3)
plot3(X_test(:,1)+R_Earth(1), X_test(:,2)+R_Earth(2), X_test(:,3)+R_Earth(3), 'm:', 'LineWidth', 2)
plot3(R_Earth(1), R_Earth(2), R_Earth(3),'b.','MarkerSize',40);
plot3(R_E_t(1,:), R_E_t(2,:), R_E_t(3,:), 'b-.', 'LineWidth', 2)
plot3(R_E_t_2(1,:), R_E_t_2(2,:), R_E_t_2(3,:), 'k:')
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]")
legend("Sun", "SC", "SC Truth Orbit", "SC Modeled Orbit", "Earth", "Earth Orbit")
view([30 35])

%% Part 1: DSN Measurement Processing
data2a = readmatrix("..\Data\Project2a_Obs.txt");
stations = readObsData(stations, data2a);






