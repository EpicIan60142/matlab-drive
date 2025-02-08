%% ASEN 5044 Project 1 Progress Check 1 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
    % Physical constants
muEarth = 398600; % km^3/s^2
rEarth = 6378; % km
omegaEarth = (2*pi/86400); % rad/s

    % Ground stations
numStations = 12; % Number of ground stations tracking the satellite
x0_s = zeros(2, numStations);
theta0_s = zeros(1, numStations);
for k = 1:numStations
    theta0_s(k) = (k-1)*(pi/6);
    x0_s(:,k) = [rEarth*cos(theta0_s(k)); rEarth*sin(theta0_s(k))];
end

    % Satellite
r0 = 6678; % km
X0 = r0; % km
n = sqrt(muEarth/(r0^3)); % rad/s
YDot0 = r0*n; % km/s

x0 = [X0; 0; 0; YDot0];
xPerturb0 = [0; 0.075; 0; -0.021];
x0 = x0 + xPerturb0;

dT = 10; % sec

T = sqrt(r0^3/muEarth); % Orbit period in seconds

%% Part 2

rhoNom = @(x_s, y_s) sqrt((X0 - x_s)^2 + y_s^2);

Abar = [
            0                   1   0               0
            2*muEarth/(X0^3)    0   0               0
            0                   0   0               1
            0                   0   -muEarth/(X0^3) 0
       ];

Bbar = [
            0   0
            1   0
            0   0
            0   1
       ];

Cbar = @(x_s, xDot_s, y_s, yDot_s) [
                                        (X0-x_s)/rhoNom(x_s, y_s)                                           0                           -y_s/rhoNom(x_s, y_s)                                                   0
                                        -y_s*(y_s*xDot_s - (X0-x_s)*(yDot0-yDot_s))/(rhoNom(x_s, y_s)^3)    (X0-x_s)/rhoNom(x_s, y_s)   (X0-x_s)*((X0-x_s)*(yDot0-yDot_s) - y_s*xDot_s)/(rhoNom(x_s, y_s)^3)    -y_s/(rhoNom(x_s, y_s)^2)
                                        y_s/(rhoNom(x_s,y_s)^2)                                             0                           (X0-x_s)/(rhoNom(x_s, y_s)^2)                                           0
                                   ];

Dbar = zeros(3,2);

Ahat = [
            Abar Bbar
            zeros(size(Bbar,2),size(Abar,2)+size(Bbar,2))
       ];

matExp = expm(Ahat*dT);

F = matExp(1:4, 1:4)
G = matExp(1:4, 5:6)
H = Cbar % Time varying - depends on x_s, xDot_s, y_s, and yDot_s!
M = Dbar

%% Part 3 - nonlinear
    % Ground station dynamics
x_stat = @(t,i) rEarth*cos(omegaEarth*t + theta0_s(i));
xDot_stat = @(t,i) -omegaEarth*rEarth*sin(omegaEarth*t + theta0_s(i));
y_stat = @(t,i) rEarth*sin(omegaEarth*t + theta0_s(i));
yDot_stat = @(t,i) omegaEarth*rEarth*cos(omegaEarth*t + theta0_s(i));
theta_stat = @(x_s, y_s) atan2(y_s, x_s);

    % Nonlinear measurements
rho_stat = @(x, x_s, y, y_s) sqrt((x-x_s)^2 + (y-y_s)^2);
rhoDot_stat = @(x, x_s, xDot, xDot_s, y, y_s, yDot, yDot_s) ...
                ((x-x_s)*(xDot-xDot_s) + (y-y_s)*(yDot-yDot_s))/rho_stat(x, x_s, y, y_s);
phi_stat = @(x, x_s, y, y_s) atan2(y-y_s, x-x_s);

    % Nonlinear dynamics simulation
tspan = 0:dT:14000;
odeOpt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t_nl, X_nl] = ode45(@(t,X)orbitEOM_unforced(t,X,muEarth), tspan, x0, odeOpt);
X_nl = X_nl';

    % Station measurements simulation
stations_nl = [];
colors = colororder('glow12');
for k = 1:length(theta0_s)
    station = struct('id', k, 't', [], 'x', [], 'xDot', [], 'y', [], 'yDot', [], 'theta', [], 'visible', [], 'rho', [], 'rhoDot', [], 'phi', []);
    station.color = colors(k,:);
    stations_nl = [stations_nl; station];
end

rho_nl = {};
rhoDot_nl = {};
phi_nl = {};
stationVis = {};
for k = 1:length(t_nl)
    x = X_nl(1,k);
    xDot = X_nl(2,k);
    y = X_nl(3,k);
    yDot = X_nl(4,k);

    for kk = 1:length(theta0_s)
        x_s = x_stat(t_nl(k),kk);
        xDot_s = xDot_stat(t_nl(k), kk);
        y_s = y_stat(t_nl(k),kk);
        yDot_s = yDot_stat(t_nl(k), kk);
        theta_s = theta_stat(x_s, y_s);

        stations_nl(kk).t = [stations_nl(kk).t; t_nl(k)];

        stations_nl(kk).x = [stations_nl(kk).x; x_s];
        stations_nl(kk).xDot = [stations_nl(kk).xDot; xDot_s];
        stations_nl(kk).y = [stations_nl(kk).y; y_s];
        stations_nl(kk).yDot = [stations_nl(kk).yDot; yDot_s];
        stations_nl(kk).theta = [stations_nl(kk).theta; theta_s];

        phi = phi_stat(x, x_s, y, y_s);

        if (phi >= (-pi/2 + theta_s)) && (phi <= (pi/2 + theta_s))
            stations_nl(kk).visible = [stations_nl(kk).visible, true];
            stations_nl(kk).rho = [stations_nl(kk).rho; rho_stat(x, x_s, y, y_s)];
            stations_nl(kk).rhoDot = [stations_nl(kk).rhoDot; rhoDot_stat(x, x_s, xDot, xDot_s, y, y_s, yDot, yDot_s)];
            stations_nl(kk).phi = [stations_nl(kk).phi; phi];
        else
            stations_nl(kk).visible = [stations_nl(kk).visible, NaN];
            stations_nl(kk).rho = [stations_nl(kk).rho; NaN];
            stations_nl(kk).rhoDot = [stations_nl(kk).rhoDot; NaN];
            stations_nl(kk).phi = [stations_nl(kk).phi; NaN];
        end

    end
end

    % Make nonlinear plots
        % Dynamics
t = tiledlayout(4,1);
title(t, "States vs. time, Nonlinear Dynamics Simulation")
nexttile;
    hold on; grid on;
    plot(t_nl, X_nl(1,:));
    xlabel("Time [sec]"); ylabel("X [km]")
nexttile;
    hold on; grid on;
    plot(t_nl, X_nl(2,:));
    xlabel("Time [sec]"); ylabel("Xdot [km/s]")
nexttile;
    hold on; grid on;
    plot(t_nl, X_nl(3,:));
    xlabel("Time [sec]"); ylabel("Y [km]")
nexttile;
    hold on; grid on;
    plot(t_nl, X_nl(4,:));
    xlabel("Time [sec]"); ylabel("Ydot [km/s]")

        % Measurements
figure; t = tiledlayout(4,1);
title(t, "Nonlinear Data Simulation")
nexttile;
    hold on; grid on;
    for k = 1:length(stations_nl)
        plot(t_nl, stations_nl(k).rho, 'x', 'Color', stations_nl(k).color);
    end
    xlabel("Time [sec]"); ylabel("rho^i [km]")
nexttile;
    hold on; grid on;
    for k = 1:length(stations_nl)
        plot(t_nl, stations_nl(k).rhoDot, 'o', 'Color', stations_nl(k).color);
    end
    xlabel("Time [sec]"); ylabel("rhoDot^i [km/s]")
nexttile;
    hold on; grid on;
    for k = 1:length(stations_nl)
        plot(t_nl, stations_nl(k).phi, 'square', 'Color', stations_nl(k).color);
    end
    xlabel("Time [sec]"); ylabel("phi^i [rad]")
nexttile;
    hold on; grid on;
    for k = 1:length(stations_nl)
        plot(t_nl, stations_nl(k).visible*stations_nl(k).id, '^', 'Color', stations_nl(k).color);
    end
    xlabel("Time [sec]"); ylabel("Visible Station ID")

%% Part 3 - Linear
t_l = t_nl;

    % Nominal trajectory (circular orbit with radius r0)
xNom_l = zeros(size(Abar,1), length(t_l));
for k = 1:length(t_l)
    xNom_l(:,k) = [
                    r0*cos(n*t_l(k));
                    -r0*n*sin(n*t_l(k));
                    r0*sin(n*t_l(k));
                    r0*n*cos(n*t_l(k))
                ];
end

    % Propagate perturbation dynamics
xPerturb_l = zeros(size(Abar,1), length(t_l));
for k = 1:length(t_l)
    xPerturb_l(:,k) = F^(k-1)*xPerturb0;
end

    % Recover total state
X_l = xNom_l + xPerturb_l;

    % Get measurements
stations_l = [];
colors = colororder('glow12');
for k = 1:length(theta0_s)
    station = struct('id', k, 't', [], 'x', [], 'xDot', [], 'y', [], 'yDot', [], 'theta', [], 'visible', [], 'rho', [], 'rhoDot', [], 'phi', []);
    station.color = colors(k,:);
    stations_l = [stations_l; station];
end

for k = 1:length(t_l)
    x = X_l(1,k);
    xDot = X_l(2,k);
    y = X_l(3,k);
    yDot = X_l(4,k);

    for kk = 1:length(theta0_s)
        x_s = x_stat(t_l(k),kk);
        xDot_s = xDot_stat(t_l(k), kk);
        y_s = y_stat(t_l(k),kk);
        yDot_s = yDot_stat(t_l(k), kk);
        theta_s = theta_stat(x_s, y_s);

        stations_l(kk).t = [stations_l(kk).t; t_l(k)];

        stations_l(kk).x = [stations_l(kk).x; x_s];
        stations_l(kk).xDot = [stations_l(kk).xDot; xDot_s];
        stations_l(kk).y = [stations_l(kk).y; y_s];
        stations_l(kk).yDot = [stations_l(kk).yDot; yDot_s];
        stations_l(kk).theta = [stations_l(kk).theta; theta_s];

        phi = phi_stat(x, x_s, y, y_s);

        if phi >= (-pi/2 + theta_s) && phi <= (pi/2 + theta_s)
            stations_l(kk).visible = [stations_l(kk).visible, true];
            stations_l(kk).rho = [stations_l(kk).rho; rho_stat(x, x_s, y, y_s)];
            stations_l(kk).rhoDot = [stations_l(kk).rhoDot; rhoDot_stat(x, x_s, xDot, xDot_s, y, y_s, yDot, yDot_s)];
            stations_l(kk).phi = [stations_l(kk).phi; phi];
        else
            stations_l(kk).visible = [stations_l(kk).visible, NaN];
            stations_l(kk).rho = [stations_l(kk).rho; NaN];
            stations_l(kk).rhoDot = [stations_l(kk).rhoDot; NaN];
            stations_l(kk).phi = [stations_l(kk).phi; NaN];
        end

    end
end

    % Make plots
        % Dynamics
figure; t = tiledlayout(4,1);
title(t, "States vs. time, Approximate Linear Dynamics Simulation")
nexttile;
    hold on; grid on;
    plot(t_l, X_l(1,:));
    xlabel("Time [sec]"); ylabel("X [km]")
nexttile;
    hold on; grid on;
    plot(t_l, X_l(2,:));
    xlabel("Time [sec]"); ylabel("Xdot [km/s]")
nexttile;
    hold on; grid on;
    plot(t_l, X_l(3,:));
    xlabel("Time [sec]"); ylabel("Y [km]")
nexttile;
    hold on; grid on;
    plot(t_l, X_l(4,:));
    xlabel("Time [sec]"); ylabel("Ydot [km/s]")

    % Perturbations
figure; t = tiledlayout(4,1);
title(t, "Approximate Linearized State Perturbations vs. time")
nexttile;
    hold on; grid on;
    plot(t_l, xPerturb_l(1,:));
    xlabel("Time [sec]"); ylabel("\deltaX [km]")
nexttile;
    hold on; grid on;
    plot(t_l, xPerturb_l(2,:));
    xlabel("Time [sec]"); ylabel("\deltaXdot [km/s]")
nexttile;
    hold on; grid on;
    plot(t_l, xPerturb_l(3,:));
    xlabel("Time [sec]"); ylabel("\deltaY [km]")
nexttile;
    hold on; grid on;
    plot(t_l, xPerturb_l(4,:));
    xlabel("Time [sec]"); ylabel("\deltaYdot [km/s]")

        % Measurements
figure; t = tiledlayout(4,1);
title(t, "Approximate Linear Data Simulation")
nexttile;
    hold on; grid on;
    for k = 1:length(stations_l)
        plot(t_l, stations_l(k).rho, 'x', 'Color', stations_l(k).color);
    end
    xlabel("Time [sec]"); ylabel("rho^i [km]")
nexttile;
    hold on; grid on;
    for k = 1:length(stations_l)
        plot(t_l, stations_l(k).rhoDot, 'o', 'Color', stations_l(k).color);
    end
    xlabel("Time [sec]"); ylabel("rhoDot^i [km]")
nexttile;
    hold on; grid on;
    for k = 1:length(stations_l)
        plot(t_l, stations_l(k).phi, 'square', 'Color', stations_l(k).color);
    end
    xlabel("Time [sec]"); ylabel("phi^i [rad]")
nexttile;
    hold on; grid on;
    for k = 1:length(stations_l)
        plot(t_l, stations_l(k).visible*stations_l(k).id, '^', 'Color', stations_l(k).color);
    end
    xlabel("Time [sec]"); ylabel("Visible Station ID")

