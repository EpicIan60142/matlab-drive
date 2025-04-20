%% ASEN 6020 Final Project Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Logistics
addpath(genpath(".\CubeSats"));
addpath(genpath(".\RaceCourse"));

%% Setup Problem Parameters
    % Define ring semi-major and minor axis ranges
semiMaj = [1, 5]; % m
semiMin = [1, 5]; % m

    % Define inter-ring parameter ranges
interRingDist = [100, 150]; % m
azimuthAngle = deg2rad([-90, 90]); % deg -> rad
elevationAngle = deg2rad([-90, 90]); % deg -> rad

    % Define number of rings range
numRings = [15, 25];

    % Define course origin orbit
T = 180*60; % 180 minute circular orbit
n = 2*pi/T;

    % Define course parameters structure
courseParams = struct("semiMaj", semiMaj, "semiMin", semiMin, "dist", interRingDist, ...
                      "azAng", azimuthAngle, "elAng", elevationAngle, "numRings", numRings, ...
                      "n", n);

    % ODE45 options
opt = odeset('AbsTol', 1e-13, 'RelTol', 1e-13);

%% Generate race course
    % Set rng seed for testing
if true
    seed = 2*69420; % cool options: 0, 2, 3, 4
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

    % Plot race course
figure(420)
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

courseCenter = scatter3(0, 0, 0, 20, 'k', 'filled', 'h');

xlabel("X [m]"); ylabel("Y [m]"); zlabel("Z [m]"); cBar = colorbar;
cBar.Label.String = "Ring Number"; colormap("cool");

%% Generate CubeSats
    % Max thrust based on https://www.cubesatshop.com/product/cluster-ifm-nano-thruster-smallsats/ - 10-500 uN
numSats = 4;
names = ["Eeny", "Meeny", "Miny", "Mo"]; % Names :P Mo has attitude control system problems, and can only use x-y-z thrusters!
markers = ["o", "^", "square", "diamond"]; % Markers for plotting
colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"]; % Colors for plotting
thrusts = [100, 250, 500, 500]*(10^-3); % Maximum thrust values
thrustConfigs = [false, false, false, true]; % Whether thrust is split amongst x-y-z thrusters or a general direction

    % Make CubeSats
cubesats = [];
for k = 1:numSats
    cubesats = [cubesats; generateCubesat(thrusts(k), thrustConfigs(k), startRing, names(k), markers(k), colors(k))];
end

    % Plot starting positions
cubeAx = []; cubeLabels = [];
for k = 1:length(cubesats)
    cubeAx = [cubeAx, plot3(cubesats(k).X0(1), cubesats(k).X0(2), cubesats(k).X0(3), 'Color', cubesats(k).color, 'Marker', cubesats(k).marker, 'MarkerFaceColor', 'k')];
    cubeLabels = [cubeLabels, sprintf("CubeSat %s", cubesats(k).name)];
end

legend([cubeStart, cubeEnd, courseCenter, cubeAx], ["CubeSat 3\sigma Starting Sample Space", "CubeSat End Target Space", "Race Course Origin", cubeLabels], 'location', 'best');

    % Show off course
for k = 0:360
    view(-30 + k, 35);
    drawnow
end

% Test dynamics
X0 = [cubesats(1).X0; 1e-3*ones(6,1)];
[t,XTest] = ode45(@(t,X)CHWEOM(t,X,cubesats(1),courseParams), [0 T], X0, opt);

% figure
% hold on; grid on;
% plot3(XTest(:,1), XTest(:,2), XTest(:,3));
% view([30 35])

% return;

%% Run the race course!
for k = 1:length(cubesats)
    fprintf("\nCubesat %s is starting the race course!\n", cubesats(k).name)
    for kk = 1:length(rings)-1 % Intermediate ring problem
            % Solve trajectory
        [t,X] = solveTrajectory_int(cubesats(k), rings(kk), courseParams, opt);
        cubesats(k).X = [cubesats(k).X; X];
        cubesats(k).t = [cubesats(k).t; t];

            % Update t0 and X0
        cubesats(k).X0 = X(end,1:6)';
        cubesats(k).t0 = t(end);

            % Plot segment
        figure
        title(sprintf("Cubesat %s trajectory segment: Ring %.0f to %.0f", cubesats(k).name, kk-1, kk));
        hold on; grid on;
        if kk == 1
            plotRing(startRing, 'g-');
        else
            plotRing(rings(kk).params.lastRing, 'g-');
        end
        plot3(X(1,1), X(1,2), X(1,3), 'k', 'Marker', cubesats(k).marker);
        plotRing(rings(kk), 'r-');
        plot3(X(:,1), X(:,2), X(:,3), 'm--');
        plot3(X(:,7), X(:,8), X(:,9), 'k--');
        view([30 35]);

            % Report progress
        fprintf("\n\tCubesat %s passed through ring %.0f in %.3f sec!", cubesats(k).name, kk, t(end) - t(1));
    end
        % Add trajectory to race plot
    figure(420)
    % title(sprintf("Cubesat %s Race Course Trajectory", cubesats(k).name));
    plot3(cubesats(k).X(:,1), cubesats(k).X(:,2), cubesats(k).X(:,3), '-', 'Color', cubesats(k).color, 'DisplayName', sprintf("Cubesat %s trajectory", cubesats(k).name));
end



