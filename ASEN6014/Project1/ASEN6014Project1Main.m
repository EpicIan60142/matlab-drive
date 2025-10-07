%% ASEN 6014 Priject 1 Main Script
% By: Ian Faber, 10/07/2025

%% Housekeeping
clc; clear; close all;

addpath("Plotting/")
addpath("RaceCourse/")
addpath("../utilities/")

%% Setup
    % Planetary constants structure
pConst.mu = 398600.4415; % km^3/m^2

    % Race course settings
        % Ring geometry
semiMaj = [1, 5]/1000; % m->km, range of ring semi-major axes
semiMin = [1, 5]/1000; % m->km, range of ring semi-minor axes

        % Course geometry
interRingDist = [300, 350]/1000; % m, range of ring center separation distances
azimuthAngle = deg2rad([-90, 90]); % deg -> rad, relative azimuth angle range between rings
elevationAngle = deg2rad([-90, 90]); % deg -> rad, relative elevation angle range between rings
numRings = 13; % Num rings - 2, need to add initial and final rings

        % Course parameters structure
courseParams = struct("semiMaj", semiMaj, "semiMin", semiMin, "dist", interRingDist, ...
                      "azAng", azimuthAngle, "elAng", elevationAngle, "numRings", numRings);

        % Course origin orbit settings
T = 180*60; % 180 minute period
oe_c.a = (pConst.mu*(T/(2*pi))^2)^(1/3);
oe_c.e = 0; 
oe_c.i = 0;
oe_c.RAAN = 0;
oe_c.argPeri = 0;
oe_c.f = 0;

    % ode45 settings
opt = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);

%% Task 1: Generate random race course
    % Determine if a random or seeded course is made
randomCourse = true;

if ~randomCourse
    seed = 3;
    rng(seed);
else
    rng("shuffle")
end

    % Generate race course
        % Intermediate rings
rings = generateRaceCourse(courseParams);

        % Starting ring
startRing = generateRing(5*max(semiMaj), 5*max(semiMin), 0, 0, -max(interRingDist), rings(1));
rings(1).params.lastRing = startRing;
rings = [startRing; rings];

        % End ring
endRing = generateRing(min(semiMaj), min(semiMin), 0, 0, min(interRingDist), rings(end));
rings = [rings; endRing];

    % Calculate race course origin and move center to be 0,0
centers = [];
for k = 1:length(rings)
    centers = [centers, rings.center];
end

newOrigin = mean(centers, 2);

for k = 1:length(rings)
    rings(k).center = rings(k).center - newOrigin;
    rings(k).params.lastRing.center = rings(k).params.lastRing.center - newOrigin;
end

    % Plot race course 
figNum = 1; 
titleText = sprintf("Generated Race Course");
xLabel = sprintf("Radial [km]"); yLabel = sprintf("Along-Track [km]"); zLabel = sprintf("Cross-Track [km]");
courseFig = plotCourse(rings, figNum, titleText, xLabel, yLabel, zLabel);

for k = 0:180
    view(-30 + 2*k, 35);
    drawnow;
end

    % Plot example ring
figNum = 2;
titleText = sprintf("Example Ring");
xLabel = sprintf("Radial [m]"); yLabel = sprintf("Along-Track [m]"); zLabel = sprintf("Cross-Track [m]");
exRingFig = plotExampleRing(rings, figNum, titleText, xLabel, yLabel, zLabel);
for k = 0:180
    view(-30 + 2*k, 35);
    drawnow;
end

 
