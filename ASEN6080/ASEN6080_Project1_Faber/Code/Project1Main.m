%% ASEN 6080 Project 1 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath(genpath("..\..\Utilities"))

earthConst = getPlanetConst();
scConst = getSCConst();
stations = makeStations();

%% Load and parse provided data
data = readmatrix("..\Data\project.txt");

stations = readData(stations, data);

%% Filter setup
stationX0 = [];
for k = 1:length(stations)
    stationX0 = [stationX0; stations(k).X0];
end

X0 = [scConst.X0_cart; earthConst.mu; earthConst.J2; scConst.Cd; stationX0];
P0 = diag([
            1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e20, 1e6, 1e6, ...
            1e-10, 1e-10, 1e-10, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6
          ]);

opt = odeset('AbsTol',1e-12,'RelTol',1e-12);

% %% Test dynamics
% [t_test, X_test] = ode45(@(t,X)orbitEOM_MuJ2Drag(t,X,earthConst,scConst), 0:20:18340, X0, opt);
% 
% Ri = earthConst.Ri; theta0 = 0;
% [earthX_s, earthY_s, earthZ_s] = sphere(20);
% earthX_s = Ri*earthX_s;
% earthY_s = Ri*earthY_s;
% earthX = earthX_s*cos(theta0) - earthY_s*sin(theta0);
% earthY = earthX_s*sin(theta0) + earthY_s*cos(theta0);
% earthZ = Ri*earthZ_s;
% I = imread('2k_earth_daymap.jpg'); % Image from https://www.solarsystemscope.com/textures/
% 
% figure; 
% hold on; grid on; axis equal;
% plot3(X_test(:,1), X_test(:,2), X_test(:,3), 'c.')
% plot3(X_test(:,10), X_test(:,11), X_test(:,12), 'b.')
% plot3(X_test(:,13), X_test(:,14), X_test(:,15), 'r.')
% plot3(X_test(:,16), X_test(:,17), X_test(:,18), 'k.')
% earth = surf(earthX, earthY, earthZ, 'FaceAlpha', 0.5);
% set(earth,'FaceColor','texturemap','cdata',I,'edgecolor','none');
% legend("Spacecraft Trajectory", "Station 101", "Station 337", "Station 394", 'Location','bestoutside');
% view([30 35])
% 
% return;

%% Run Batch
batchRun = runBatch(X0, zeros(size(X0)), P0, earthConst, scConst, stations, data(:,1), [], opt, 3);


%% Run LKF
LKFRun = runLKF(X0, zeros(size(X0)), P0, earthConst, scConst, stations, 10);



