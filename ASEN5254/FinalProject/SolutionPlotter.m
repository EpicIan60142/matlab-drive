%% ASEN 5254 Final Project Solution Plotter
% By: Ian Faber, 12/03/2025

%% Housekeeping
clc; clear; close all;

%% Setup
% Choose dataset
load("Run_05-May-2025_215651_Scenario2.mat");

% Read in solutions
for k = 1:4
    if k == 1
        solution = readmatrix("Eeny.txt");
    elseif k == 2
        solution = readmatrix("Meeny.txt");
    elseif k == 3
        solution = readmatrix("Miny.txt");
    else
        solution = readmatrix("Mo.txt");
    end

    cubesats(k).X_PVT = cubesats(k).X;
    cubesats(k).t_PVT = cubesats(k).t;
    cubesats(k).u_PVT = cubesats(k).u;
    cubesats(k).X_SST = solution(:,1:6);
    cubesats(k).t_SST = cumsum(solution(:,end));
    cubesats(k).u_SST = solution(:,7:9);
    
    cubesats(k).X = cubesats(k).X_SST;
    cubesats(k).t = cubesats(k).t_SST;
    cubesats(k).u = cubesats(k).u_SST;
end

%% Plot race course
fig = plotCourse(startRing, rings, endRing, cubesats, 1, "Race Course", "Radial [m]", "Along-Track [m]", "Cross-Track [m]", true);

%% Plot trajectory differences
for k = 1:4
    plotFullTrajectory(cubesat(k), rings, k+1, "Cubesat Trajectory");
end



