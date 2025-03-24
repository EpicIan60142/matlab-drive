%% ASEN 6020 Final Project Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup Problem Parameters
    % Define ring semi-major and minor axis ranges
semiMaj = [1, 5]; % m
semiMin = [1, 5]; % m

    % Define inter-ring parameter ranges
interRingDist = [15, 50]; % m
azimuthAngle = deg2rad([-90, 90]); % deg -> rad
elevationAngle = deg2rad([-90, 90]); % deg -> rad

    % Define number of rings range
numRings = [15, 25];

    % Define ring parameters structure
courseParams = struct("semiMaj", semiMaj, "semiMin", semiMin, "dist", interRingDist, ...
                      "azAng", azimuthAngle, "elAng", elevationAngle, "numRings", numRings);

%% Generate race course
    % Set rng seed for testing
if true
    seed = 2*69420;
    rng(seed);
end

    % Make course
rings = generateRaceCourse(courseParams);

    % Define CubeSat starting ring as the maximum distance behind the start
    % ring at a 3 sigma covariance of 5x the largest semi-major/minor axis
startRing = generateRing(50*max(semiMaj), 50*max(semiMin), 0, 0, -max(interRingDist), rings(1));

    % Define CubeSat ending ring as the minimum distance after the last 
    % ring at a covariance corresponding to the smallest possible ring
endRing = generateRing(min(semiMaj),min(semiMin),0,0,min(interRingDist),rings(end));
rings = [rings; endRing];

    % Plot course
figure
hold on; grid on; axis equal
title("Generated Race Course")
for k = 1:length(rings)-1
    scatter3(rings(k).center(1), rings(k).center(2), rings(k).center(3), 10, k, 'filled')
    quiver3(rings(k).center(1), rings(k).center(2), rings(k).center(3), rings(k).normal(1), rings(k).normal(2), rings(k).normal(3), 10, 'filled', 'k-')
    plotRing(rings(k), 'k-');
end
cubeStart = plotRing(startRing, 'g-');
quiver3(startRing.center(1), startRing.center(2), startRing.center(3), startRing.normal(1), startRing.normal(2), startRing.normal(3), 10, 'filled', 'k-')
cubeEnd = plotRing(endRing, 'r-');

courseCenter = scatter3(0, 0, 0, 10, 'k', 'filled');

xlabel("X [m]"); ylabel("Y [m]"); zlabel("Z [m]"); cBar = colorbar;
cBar.Label.String = "Ring Number"; colormap("cool");

legend([cubeStart, cubeEnd, courseCenter], ["CubeSat 3\sigma Starting Sample Space", "CubeSat End Target Space", "Race Course Origin"], 'location', 'best');

for k = 0:360
    view(-30 + k, 35);
    drawnow
end
