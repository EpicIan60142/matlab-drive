%% Ian Faber - ASEN 2004 Rocket Equation Model
%--------------------------------------------------------------------------
%   By: Ian Faber
%   SID: 108577813
%   Started: 4/3/22, 6:25 PM
%   Finished: 4/11/22, 9:55 PM
%   
%   Runs a bottle rocket simulation subject to a set of initial parameters
%   according to delta-V calculated from the rocket equation with ode45 in 
%   3 dimensions. The simulation includes wind, but only to control the
%   heading of the rocket (no side forces from wind). A Monte Carlo
%   simulation is alo run to analyze uncertainty in the rocket's landing
%   coordinates.
%
%--------------------------------------------------------------------------

%% Housekeeping
clc; clear; close all;

% Choose LA trial to reference
rocketTrial = 1;

%% Setup - Single flight

const = getConst(rocketTrial, 0, 0, 0, 0, 1);

% Calculate initial x, y, and z velocities
vx0 = const.deltaV*cosd(const.thetaInit);
vy0 = 0;
vz0 = const.deltaV*sind(const.thetaInit);

% Format the initial conditions vector, and by extension the variables to
% integrate
X0 = [const.xInit; const.yInit; const.zInit; vx0; vy0; vz0];

% Define events worthy of stopping integration, i.e. hitting the ground
options = odeset('Events', @phase);

%% Simulation - Single flight

% Integrate! Solves for the trajectory of the rocket by integrating the
% variables in X0 over tspan according to the derivative information
% contained in rocketEOM. Also stops integration according to "options," a
% predefined set of stopping conditions
[time, state, timePhases, ~, ~] = ode45(@(t,state)rocketEOM(t,state,const, 1), const.tspan, X0, options);

% Extract intermediate variables from rocketEOM for debugging, particularly
% weight, drag, friction, and wind. Found this approach on the MATLAB
% forums.
[~,gravCell, dragCell, fricCell, windCell] = cellfun(@(t,state)rocketEOM(t,state.',const, 1), num2cell(time), num2cell(state,2), 'uni', 0);

%Allocate space for intermediate variables
gravity = zeros(length(time),1);
drag = zeros(length(time),1);
friction = zeros(length(time),1);
wind = zeros(length(time),1);

% Extract intermediate variables from their cells
for i = 1:length(time)
    gravity(i) = norm(gravCell{i});
    drag(i) = norm(dragCell{i});
    friction(i) = norm(fricCell{i});
    wind(i) = norm(windCell{i});
end

%% Setup - Monte Carlo simulation

% Run 100 cases
nSims = 100;

% Generate a new contant structure with vectors for each simulation
monteConst = getConst(rocketTrial, 0.0005, 0.0005, 1, 1, nSims);

% Calculate initial velocities
monteVx0 = monteConst.deltaV.*cosd(monteConst.thetaInit);
monteVy0 = zeros(nSims, 1);
monteVz0 = monteConst.deltaV.*sind(monteConst.thetaInit);

%% Simulation - Monte Carlo

% Preallocate structures and arrays for speed
monteX = struct([]);
monteY = struct([]);
monteZ = struct([]);
monteRange = zeros(nSims, 1);
monteCrossRange = zeros(nSims, 1);
monteHeight = zeros(nSims, 1);

% Run simulations with randomized values for various parameters
for k = 1:nSims
    monteX0 = [monteConst.xInit, monteConst.yInit, monteConst.zInit, monteVx0(k), monteVy0(k), monteVz0(k)];
    [monteTime, monteState, monteTimePhases, ~, ~] = ode45(@(t,state)rocketEOM(t,state,monteConst, k), monteConst.tspan, monteX0, options);

    monteX{k} = monteState(:,1);
    monteY{k} = monteState(:,2);
    monteZ{k} = monteState(:,3);
    
    [~, maxR] = max(abs(monteX{k}));
    [~, maxCR] = max(abs(monteY{k}));
    [~, maxH] = max(abs(monteZ{k}));

    monteRange(k) = monteX{k}(maxR);
    monteCrossRange(k) = monteY{k}(maxCR);
    monteHeight(k) = monteZ{k}(maxH);
end

%% Extraction - Single Flight

% Extract variables of interest
rocketX = state(:,1);
rocketY = state(:,2);
rocketZ = state(:,3);
rocketVx = state(:,4);
rocketVy = state(:,5);
rocketVz = state(:,6);

velocityNorms = sqrt(rocketVx.^2 + rocketVy.^2 + rocketVz.^2);

% Find maximum values of interest
maxRange = max(rocketX);
maxCrossRange = max(abs(rocketY));
distance = norm([maxRange, maxCrossRange]);
maxHeight = max(rocketZ);
maxVx = max(rocketVx);
maxVy = max(abs(rocketVy));
maxVz = max(rocketVz);

alpha = atand(maxCrossRange/maxRange);

northLine = @(x) (sind(const.thetaAim)/cosd(const.thetaAim))*x;
eastLine = @(x) (sind(const.thetaAim + 90)/cosd(const.thetaAim + 90))*x;

[~, thetaAloft, thetaGround] = analyzeWind(const, const.hAloft, 1);

windLine = @(x) (sind(const.thetaAim - thetaAloft)/cosd(const.thetaAim - thetaAloft))*x;

%% Plotting

% Plot the trajectory and variables of interest for the bottle rocket's
% flight!
f = figure();
%f.Position = [100 100 740 740];

% Single Flight Trajectory
%subplot(2,2,1)
hold on;
title("Bottle Rocket Full Trajectory");
color_line3d(time, rocketX, rocketY, rocketZ);

% Create arrows showing cardinal and wind directions
mArrow3([0, northLine(0), 0], [80, northLine(80), 0], 'color', 'r', 'stemWidth', 0.25, 'tipWidth', 1);
mArrow3([0, northLine(0), 0], [80, eastLine(80), 0], 'color', 'g', 'stemWidth', 0.25, 'tipWidth', 1);
mArrow3([0, northLine(0), 0], [-80, northLine(-80), 0], 'color', 'b', 'stemWidth', 0.25, 'tipWidth', 1);
mArrow3([0, northLine(0), 0], [-80, eastLine(-80), 0], 'color', 'm', 'stemWidth', 0.25, 'tipWidth', 1);
mArrow3([-80, windLine(-80), 0], [80, windLine(80), 0], 'color', 'c', 'stemWidth', 0.25, 'tipWidth', 1);

% Plot x and y axes
plot3(-80:1:80, zeros(161,1), zeros(161,1), 'k--');
plot3(zeros(161,1), -80:1:80, zeros(161,1), 'k-.');

xlim([-90, 90]);
ylim([-90, 90]);
zlim([0, 30]);
view([0 90]);
%view([30 35]);
xlabel("Range (m)");
ylabel("Crossrange (m)");
zlabel("Height (m)");

windLabel = sprintf("Wind aloft: from %s", const.windDirAloft);

legend("Trajectory", "North", "East", "South", "West", windLabel, "x Axis", "y axis")

hold off;

%{
% Drag
subplot(2,2,2)
hold on;
title("Bottle Rocket Drag Force");
plot(time, drag);
xlabel("Time (sec)");
ylabel("Drag (N)");
hold off;

% Wind
subplot(2,2,3)
hold on;
title("Bottle Rocket Wind");
plot(time, wind);
xlabel("Time (sec)");
ylabel("Windspeed (m/s)");
hold off;

% Friction
endTime = 40;
subplot(2,2,4)
hold on;
title("Bottle Rocket Friction from Stand");
plot(time(1:endTime), friction(1:endTime));
xlabel("Time (sec)");
ylabel("Friction (N)");
hold off;
%}


% Monte Carlo Coordinates
g = figure();

% Plot raw landing positions
plot(monteRange,monteCrossRange,'k.','markersize',6)
axis equal;
grid on;
title("Monte Carlo Landing Coordinates")
xlabel('Range [m]');
ylabel('Crossrange [m]');
hold on;
 
% Given error ellipse code below, substituted "x" for "monteRange" and "y"
% for "monteCrossRange"

% Calculate covariance matrix
P = cov(monteRange,monteCrossRange);
mean_x = mean(monteRange);
mean_y = mean(monteCrossRange);
 
% Calculate the define the error ellipses
n=100; % Number of points around ellipse
p=0:pi/n:2*pi; % angles around a circle
 
[eigvec,eigval] = eig(P); % Compute eigen-stuff
xy_vect = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x_vect = xy_vect(:,1);
y_vect = xy_vect(:,2);
 
% Plot the error ellipses overlaid on the same figure
sigma1 = plot(1*x_vect+mean_x, 1*y_vect+mean_y, 'b');
sigma2 = plot(2*x_vect+mean_x, 2*y_vect+mean_y, 'g');
sigma3 = plot(3*x_vect+mean_x, 3*y_vect+mean_y, 'r');

legend([sigma1, sigma2, sigma3], "1\sigma", "2\sigma", "3\sigma")

% Monte Carlo Trajectories
h = figure();

title("Monte Carlo Trajectories")

% Plot all simulated trajectories
for k = 1:nSims
    hold on;
    plot3(monteX{k}, monteY{k}, monteZ{k})
end

% Overlay erorr ellipses on landing
plot(monteRange,monteCrossRange,'k.','markersize',6);
sigma1 = plot(1*x_vect+mean_x, 1*y_vect+mean_y, 'b');
sigma2 = plot(2*x_vect+mean_x, 2*y_vect+mean_y, 'g');
sigma3 = plot(3*x_vect+mean_x, 3*y_vect+mean_y, 'r');

legend([sigma1, sigma2, sigma3], "1\sigma", "2\sigma", "3\sigma")

xlabel("Range (m)");
ylabel("Crossrange (m)");
zlabel("Height (m)");
view([30 35]);

hold off;

