%% ASEN 5090 HW1 Problem 2 Main Script
% By: Ian Faber, 08/29/2025

%% Housekeeping
clc; clear; close all;

%% Setup
    % Aircraft
x0 = 250; % m
v0 = 50; % m
h0 = 100; % m

    % Ground station
x_R = 0; % m
y_R = 0; % m

%% Part a. Calculating range, range rate, and zenith angle

% Construct array of transmitter states
x_T = -x0:1:x0; % intervals of 1 m

% Calculate range, range rate, and zenith angle
R = sqrt((x_R - x_T).^2 + (y_R - h0)^2);
RDot = -((x_R - x_T)./R)*v0; % No vertical velocity
zenith = asin(x_T./R);

% Plot results
fig = figure; tl = tiledlayout(3,1); ax = [];
title(tl, "Ground Station Measurements to Aircraft");
nt = nexttile; ax(end+1) = nt;
    hold on; grid on;
    title("Range vs. X_T")
    plot(x_T, R, 'b-');
    xlabel("X_T [m]"); ylabel("R [m]", 'Interpreter','latex');
nt = nexttile; ax(end+1) = nt;
    hold on; grid on;
    title("Range Rate vs. X_T")
    plot(x_T, RDot, 'r-');
    xlabel("X_T [m]"); ylabel('$\dot{R}$ [m/s]', 'Interpreter','latex');
nt = nexttile; ax(end+1) = nt;
    hold on; grid on;
    title("Zenith Angle vs. X_T")
    plot(x_T, rad2deg(zenith), 'k-');
    xlabel("X_T [m]"); ylabel('$\xi$ [deg]', 'Interpreter','latex');
linkaxes(ax, 'x');

%% Part b-c. Sensitivity analysis

% Test variations on x and y
dx = -10:10; % horizontal deviation bounds
dy = -10:10; % vertical deviation bounds

[dX, dY] = meshgrid(dx, dy);

% Calculate dR, dRDot, and dZenith at designated test points on the
% aircraft trajectory
x_test = x_T;
dR = []; dRDot = []; dZenith = [];

for k = 1:length(x_test)
    R = sqrt((x_R - x_test(k)).^2 + (y_R - h0)^2);

    dR(:,:,k) = -((x_R - x_test(k))/R)*dX - ((y_R - h0)/R)*dY;
    dRDot(:,:,k) = (v0/R)*dX + (((x_R - x_test(k))*v0)/(R^2))*dR(:,:,k);
    dZenith(:,:,k) = (R*dX - x_test(k)*dR(:,:,k))/((R^2)*sqrt(1-(x_test(k)/R)^2));

end

% Plot Range results
numLevels = 20; % Number of contours to plot
x_plot = floor(linspace(-x0, x0, 9));
midpoint_dx = find(dx == 0);
midpoint_dy = find(dy == 0);
bounds = 2; % How many deviations outside of 0 to plot sensitivity vs. x_T

idx = [];
for k = 1:length(x_plot)
    num = find(x_plot(k) == x_test);

    if ~isempty(num)
        idx(end+1) = num;
    end
end

fig = figure; tl = tiledlayout('flow'); ax = [];
for k = idx
    mindR = min(dR(:,:,k), [], 'all');
    maxdR = max(dR(:,:,k), [], 'all');
    levels = linspace(mindR, maxdR, numLevels);
    nt = nexttile; ax(end+1) = nt;
        hold on; grid on;
        title(sprintf("Range Sensitivity vs. aircraft position at x_T = %.3f m", x_test(k)), 'FontSize', 9)
        contourf(dX, dY, dR(:,:,k), levels)
        c = colorbar; c.Label.String = "\deltaR [m]"; c.Location = "southoutside";
        xlabel("\deltax [m]"); ylabel("\deltay [m]")
end
linkaxes(ax, 'x', 'y');
fig.Position = [10 10 1200 1000];

% Plot Range Rate results
fig = figure; tl = tiledlayout('flow'); ax = [];
for k = idx
    mindRDot = min(dRDot(:,:,k), [], 'all');
    maxdRDot = max(dRDot(:,:,k), [], 'all');
    levels = linspace(mindRDot, maxdRDot, numLevels);
    nt = nexttile; ax(end+1) = nt;
        hold on; grid on;
        title(sprintf("Range Rate Sensitivity vs. aircraft position at x_T = %.3f m", x_test(k)), 'FontSize', 8.5)
        contourf(dX, dY, dRDot(:,:,k), levels)
        c = colorbar; c.Label.Interpreter = 'latex'; c.Label.String = "$\delta \dot{R}$ [m/s]";  c.Location = "southoutside";
        xlabel("\deltax [m]"); ylabel("\deltay [m]")
end
linkaxes(ax, 'x', 'y');
fig.Position = [10 10 1200 1000];

% Plot Zenith Angle results
fig = figure; tl = tiledlayout('flow'); ax = [];
for k = idx
    mindZenith = min(dZenith(:,:,k), [], 'all');
    maxdZenith = max(dZenith(:,:,k), [], 'all');
    levels = linspace(mindZenith, maxdZenith, numLevels);
    nt = nexttile; ax(end+1) = nt;
        hold on; grid on;
        title(sprintf("Zenith Angle Sensitivity vs. aircraft position at x_T = %.3f m", x_test(k)), 'FontSize', 8.5)
        contourf(dX, dY, rad2deg(dZenith(:,:,k)), rad2deg(levels))
        c = colorbar; c.Label.String = "\delta\xi [deg]"; c.Location = "southoutside";
        xlabel("\deltax [m]"); ylabel("\deltay [m]")
end
linkaxes(ax, 'x', 'y');
fig.Position = [10 10 1200 1000];

figure; plots = []; labels = [];
hold on; grid on;
title("Range Sensitivity vs. x_T")
for k = -bounds:bounds
    plots = [plots; plot(x_test, squeeze(dR(midpoint_dx+k,midpoint_dy+k,:)))];
    labels = [labels; sprintf("\\deltax = %.0f m, \\deltay = %.0f m", dx(midpoint_dx+k), dy(midpoint_dy+k))];
end
xlabel("x_T [m]"); ylabel("\deltaR [m]")
legend(plots, labels, 'location', 'eastoutside')

figure; plots = []; labels = [];
hold on; grid on;
title("Range Rate Sensitivity vs. x_T")
for k = -bounds:bounds
    plots = [plots; plot(x_test, squeeze(dRDot(midpoint_dx+k,midpoint_dy+k,:)))];
    labels = [labels; sprintf("\\deltax = %.0f m, \\deltay = %.0f m", dx(midpoint_dx+k), dy(midpoint_dy+k))];
end
xlabel("x_T [m]"); ylabel("$\delta\dot{R}$ [m/s]", 'interpreter', 'latex')
legend(plots, labels, 'location', 'eastoutside')

figure; plots = []; labels = [];
hold on; grid on;
title("Zenith Angle Sensitivity vs. x_T")
for k = -bounds:bounds
    plots = [plots; plot(x_test, squeeze(rad2deg(dZenith(midpoint_dx+k,midpoint_dy+k,:))))];
    labels = [labels; sprintf("\\deltax = %.0f m, \\deltay = %.0f m", dx(midpoint_dx+k), dy(midpoint_dy+k))];
end
xlabel("x_T [m]"); ylabel("\delta\xi [deg]")
legend(plots, labels, 'location', 'eastoutside')


