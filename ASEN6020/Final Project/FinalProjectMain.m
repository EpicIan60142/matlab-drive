%% ASEN 6020 Final Project Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup Problem Parameters
    % Define ring semi-major and minor axis ranges
semiMaj = [1, 3]; % m
semiMin = [1, 3]; % m

    % Define inter-ring parameter ranges
interRingDist = [25, 50]; % m
azimuthAngle = deg2rad([-60, 60]); % deg -> rad
elevationAngle = deg2rad([-60 60]); % deg -> rad

    % Define number of rings range
numRings = [15, 25];

    % Define ring parameters structure
courseParams = struct("semiMaj", semiMaj, "semiMin", semiMin, "dist", interRingDist, ...
                      "azAng", azimuthAngle, "elAng", elevationAngle, "numRings", numRings);

%% Generate race course
seed = 69420;
rng(seed); % Set rng seed for testing

rings = generateRaceCourse(courseParams);

figure
hold on; grid on; axis equal
title("Generated Race Course")
for k = 1:length(rings)
    scatter3(rings(k).center(1), rings(k).center(2), rings(k).center(3),15,k,'filled')
    quiver3(rings(k).center(1), rings(k).center(2), rings(k).center(3), rings(k).normal(1), rings(k).normal(2), rings(k).normal(3), 10, 'filled', 'k-')
end
xlabel("X [m]"); ylabel("Y [m]"); zlabel("Z [m]"); cBar = colorbar;
view([-30 35]); cBar.Label.String = "Ring Number";
