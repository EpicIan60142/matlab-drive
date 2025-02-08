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

xNom0 = [X0; 0; 0; YDot0];
xPerturb0 = [0; 0.075; 0; -0.021];
x0 = xNom0 + xPerturb0;

dT = 10; % sec

T = 2*pi*sqrt(r0^3/muEarth); % Orbit period in seconds

%% Part 2

rhoNom = @(x, x_s, y, y_s) sqrt((x-x_s)^2 + (y-y_s)^2);

Atilde = @(x, y) ...
    [
        0                                        1   0                                               0
        muEarth*(2*x^2-y^2)/((x^2+y^2)^(5/2))    0   3*(muEarth*x*y)/((x^2+y^2)^(5/2))               0
        0                                        0   0                                               1
        3*(muEarth*x*y)/((x^2+y^2)^(5/2))        0   muEarth*(2*y^2-x^2)/((x^2+y^2)^(5/2))           0
    ];

Btilde = [
            0   0
            1   0
            0   0
            0   1
         ];

Ctilde = @(x, x_s, xDot, xDot_s, y, y_s, yDot, yDot_s) ... 
    [
         (x-x_s)/rhoNom(x, x_s, y, y_s)                                                        0                                (y-y_s)/rhoNom(x, x_s, y, y_s)                                                        0
         (y-y_s)*((y-y_s)*(xDot-xDot_s) - (x-x_s)*(yDot-yDot_s))/(rhoNom(x, x_s, y, y_s)^3)    (x-x_s)/rhoNom(x, x_s, y, y_s)   (x-x_s)*((x-x_s)*(yDot-yDot_s) - (y-y_s)*(xDot-xDot_s))/(rhoNom(x, x_s, y, y_s)^3)    (y-y_s)/(rhoNom(x, x_s, y, y_s))
         -(y-y_s)/(rhoNom(x, x_s, y, y_s)^2)                                                   0                                (x-x_s)/(rhoNom(x, x_s, y, y_s)^2)                                                    0
    ];

Dtilde = zeros(3,2);

Ahat = @(x,y) ...
    [
        Atilde(x,y) Btilde
        zeros(size(Btilde,2),size(Atilde(x,y),2)+size(Btilde,2))
    ];

matExp = @(x,y) expm(Ahat(x,y)*dT);

F = @(x,y) eye(size(Atilde(x,y))) + dT*Atilde(x,y); % Time varying - depends on x and y!
G = dT*Btilde; % Time invariant
H = Ctilde; % Time varying - depends on x, x_s, xDot, xDot_s, y, y_s, yDot, and yDot_s!
M = Dtilde; % Time invariant

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
        
        if (y-y_s) < 0 && y > 0 && x < 0 % Quadrant 2
            phi_check = phi + 2*pi;
        elseif (y-y_s) > 0 && y < 0 && x < 0 % Quadrant 3
            phi_check = phi - 2*pi;
        else
            phi_check = phi;
        end

        viewingAngle = wrapToPi(phi_check - theta_s);

        if (viewingAngle >= -pi/2) && (viewingAngle <= pi/2) %(phi_check >= (-pi/2 + theta_s)) && (phi_check <= (pi/2 + theta_s))
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
checkPoints = [0, 1350, 4010, 7250, 10970, 14000];
idxPoints = t_nl == checkPoints;
        % Dynamics
t = tiledlayout(4,1);
title(t, "States vs. time, Nonlinear Dynamics Simulation")
nexttile;
    hold on; grid on;
    nt = plot(t_nl, X_nl(1,:));
    xlabel("Time [sec]"); ylabel("X [km]")
    for k = 1:length(checkPoints)
        datatip(nt, t_nl(idxPoints(:,k)), X_nl(1,idxPoints(:,k)));
    end
nexttile;
    hold on; grid on;
    nt = plot(t_nl, X_nl(2,:));
    xlabel("Time [sec]"); ylabel("Xdot [km/s]")
    for k = 1:length(checkPoints)
        datatip(nt, t_nl(idxPoints(:,k)), X_nl(2,idxPoints(:,k)));
    end
nexttile;
    hold on; grid on;
    nt = plot(t_nl, X_nl(3,:));
    xlabel("Time [sec]"); ylabel("Y [km]")
    for k = 1:length(checkPoints)
        datatip(nt, t_nl(idxPoints(:,k)), X_nl(3,idxPoints(:,k)));
    end
nexttile;
    hold on; grid on;
    nt = plot(t_nl, X_nl(4,:));
    xlabel("Time [sec]"); ylabel("Ydot [km/s]")
    for k = 1:length(checkPoints)
        datatip(nt, t_nl(idxPoints(:,k)), X_nl(4,idxPoints(:,k)));
    end

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

    if false
        idxProblem = 290; station = 7;
        checkSituation(X_nl, stations_nl, idxProblem, station);
        
        idxProblem = 336; station = 8;
        checkSituation(X_nl, stations_nl, idxProblem, station);
        
        idxProblem = 376; station = 9;
        checkSituation(X_nl, stations_nl, idxProblem, station);
        
        idxProblem = 818; station = 6;
        checkSituation(X_nl, stations_nl, idxProblem, station);
        
        idxProblem = 863; station = 7;
        checkSituation(X_nl, stations_nl, idxProblem, station);
        
        idxProblem = 906; station = 8;
        checkSituation(X_nl, stations_nl, idxProblem, station);
    end
    
    if false
        fprintf("\n");
    
        for idxProblem = 266:274
            station = 7;
            checkSituation(X_nl, stations_nl, idxProblem, station)
        end
        
        fprintf("\n");
    
        for idxProblem = 793:813
            station = 6;
            checkSituation(X_nl, stations_nl, idxProblem, station)
        end
    
        fprintf("\n");
    
        for idxProblem = 1349:1370
            station = 5;
            checkSituation(X_nl, stations_nl, idxProblem, station)
        end
    end

%% Part 3 - Linear
t_l = t_nl;

    % Nominal trajectory (circular orbit with radius r0)
xNom_l = zeros(size(Btilde,1), length(t_l));
for k = 1:length(t_l)
    xNom_l(:,k) = [
                        r0*cos(n*t_l(k));
                        -r0*n*sin(n*t_l(k));
                        r0*sin(n*t_l(k));
                        r0*n*cos(n*t_l(k))
                  ];
end
% xNom_l = X_nl;

    % Propagate perturbation dynamics
xPerturb_l = zeros(size(Btilde,1), length(t_l));
xPerturb = xPerturb0;
xPerturb_l(:,1) = xPerturb;
for k = 2:length(t_l)
    x = xNom_l(1,k-1);
    y = xNom_l(3,k-1);
    xPerturb_l(:,k) = F(x,y)*xPerturb;
    xPerturb = xPerturb_l(:,k);
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

        rhoNom = rho_stat(xNom_l(1,k), x_s, xNom_l(3,k), y_s);
        rhoDotNom = rhoDot_stat(xNom_l(1,k), x_s, xNom_l(2,k), xDot_s, xNom_l(3,k), y_s, xNom_l(4,k), yDot_s);
        phiNom = phi_stat(xNom_l(1,k), x_s, xNom_l(3,k), y_s);
        
        meas = [rhoNom; rhoDotNom; phiNom] + H(xNom_l(1,k), x_s, xNom_l(2,k), xDot_s, xNom_l(3,k), y_s, xNom_l(4,k), yDot_s)*xPerturb_l(:,k);

        if k == 1 && kk == 1 && true
            H(xNom_l(1,k), x_s, xNom_l(2,k), xDot_s, xNom_l(3,k), y_s, xNom_l(4,k), yDot_s)
            xPerturb_l(:,k)
            H(xNom_l(1,k), x_s, xNom_l(2,k), xDot_s, xNom_l(3,k), y_s, xNom_l(4,k), yDot_s)*xPerturb_l(:,k)
            meas
        end

        rho = meas(1);
        rhoDot = meas(2);
        phi = meas(3);

        % if (y-y_s) < 0 && y > 0 && x < 0 % Quadrant 2
        %     phi_check = phi + 2*pi;
        % elseif (y-y_s) > 0 && y < 0 && x < 0 % Quadrant 3
        %     phi_check = phi - 2*pi;
        % else
        %     phi_check = phi;
        % end
        % 
        % viewingAngle = wrapToPi(phi_check - theta_s);

        if stations_nl(kk).visible(k) == true
            stations_l(kk).visible = [stations_l(kk).visible, true];
            stations_l(kk).rho = [stations_l(kk).rho; rho];
            stations_l(kk).rhoDot = [stations_l(kk).rhoDot; rhoDot];
            stations_l(kk).phi = [stations_l(kk).phi; phi];
        else
            stations_l(kk).visible = [stations_l(kk).visible, NaN];
            stations_l(kk).rho = [stations_l(kk).rho; NaN];
            stations_l(kk).rhoDot = [stations_l(kk).rhoDot; NaN];
            stations_l(kk).phi = [stations_l(kk).phi; NaN];
        end

        % if (viewingAngle >= -pi/2) && (viewingAngle <= pi/2) %(phi_check >= (-pi/2 + theta_s)) && (phi_check <= (pi/2 + theta_s))
        %     stations_l(kk).visible = [stations_l(kk).visible, true];
        %     stations_l(kk).rho = [stations_l(kk).rho; rho];
        %     stations_l(kk).rhoDot = [stations_l(kk).rhoDot; rhoDot];
        %     stations_l(kk).phi = [stations_l(kk).phi; phi];
        % else
        %     stations_l(kk).visible = [stations_l(kk).visible, NaN];
        %     stations_l(kk).rho = [stations_l(kk).rho; NaN];
        %     stations_l(kk).rhoDot = [stations_l(kk).rhoDot; NaN];
        %     stations_l(kk).phi = [stations_l(kk).phi; NaN];
        % end

    end
end

    % Make plots
checkPoints = [0, 2840, 6560, 9620, 14000];
idxPoints = t_l == checkPoints;
        % Dynamics
figure; t = tiledlayout(4,1);
title(t, "States vs. time, Approximate Linear Dynamics Simulation")
nexttile;
    hold on; grid on;
    nt = plot(t_l, X_l(1,:));
    xlabel("Time [sec]"); ylabel("X [km]")
    for k = 1:length(checkPoints)
        datatip(nt, t_l(idxPoints(:,k)), X_l(1,idxPoints(:,k)));
    end
nexttile;
    hold on; grid on;
    nt = plot(t_l, X_l(2,:));
    xlabel("Time [sec]"); ylabel("Xdot [km/s]")
    for k = 1:length(checkPoints)
        datatip(nt, t_l(idxPoints(:,k)), X_l(2,idxPoints(:,k)));
    end
nexttile;
    hold on; grid on;
    nt = plot(t_l, X_l(3,:));
    xlabel("Time [sec]"); ylabel("Y [km]")
    for k = 1:length(checkPoints)
        datatip(nt, t_l(idxPoints(:,k)), X_l(3,idxPoints(:,k)));
    end
nexttile;
    hold on; grid on;
    nt = plot(t_l, X_l(4,:));
    xlabel("Time [sec]"); ylabel("Ydot [km/s]")
    for k = 1:length(checkPoints)
        datatip(nt, t_l(idxPoints(:,k)), X_l(4,idxPoints(:,k)));
    end

    % Perturbations
checkPoints = [0, 3470, 7410, 11060, 14000];
idxPoints = t_l == checkPoints;
figure; t = tiledlayout(4,1);
title(t, "Approximate Linearized State Perturbations vs. time")
nexttile;
    hold on; grid on;
    nt = plot(t_l, xPerturb_l(1,:));
    xlabel("Time [sec]"); ylabel("\deltaX [km]")
    for k = 1:length(checkPoints)
        datatip(nt, t_l(idxPoints(:,k)), xPerturb_l(1,idxPoints(:,k)));
    end
nexttile;
    hold on; grid on;
    nt = plot(t_l, xPerturb_l(2,:));
    xlabel("Time [sec]"); ylabel("\deltaXdot [km/s]")
    for k = 1:length(checkPoints)
        datatip(nt, t_l(idxPoints(:,k)), xPerturb_l(2,idxPoints(:,k)));
    end
nexttile;
    hold on; grid on;
    nt = plot(t_l, xPerturb_l(3,:));
    xlabel("Time [sec]"); ylabel("\deltaY [km]")
    for k = 1:length(checkPoints)
        datatip(nt, t_l(idxPoints(:,k)), xPerturb_l(3,idxPoints(:,k)));
    end
nexttile;
    hold on; grid on;
    nt = plot(t_l, xPerturb_l(4,:));
    xlabel("Time [sec]"); ylabel("\deltaYdot [km/s]")
    for k = 1:length(checkPoints)
        datatip(nt, t_l(idxPoints(:,k)), xPerturb_l(4,idxPoints(:,k)));
    end

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

% for idxProblem = 918:937
%     station = 8;
%     checkSituation(X_l, stations_l, idxProblem, station);
% end

% for idxProblem = 1158:1163
%     station = 1;
%     checkSituation(X_l, stations_l, idxProblem, station);
% end
% 
% for idxProblem = 1061:1072
%     station = 11;
%     checkSituation(X_l, stations_l, idxProblem, station);
% end
% 
% for idxProblem = 1109:1117
%     station = 12;
%     checkSituation(X_l, stations_l, idxProblem, station);
% end

%% Utilities
function checkSituation(X, stations, idxProblem, station)
    idx = (idxProblem-10):(idxProblem+10);

    figure; axis equal;
    hold on; grid on;
    sat = plot(X(1,idx), X(3,idx), 'b-', 'LineWidth', 3);
    plot(X(1,:), X(3,:), 'b.', 'MarkerSize', 1);
    plot(X(1,idxProblem), X(3,idxProblem), 'k.', 'MarkerSize', 15)
    plot(X(1,idxProblem-1), X(3,idxProblem-1), 'g.', 'MarkerSize', 15)
    stat = plot(stations(station).x(idx), stations(station).y(idx), 'r-', 'LineWidth', 3);
    plot(stations(station).x(1:5:end), stations(station).y(1:5:end), 'r.', 'MarkerSize', 1)
    plot(stations(station).x(idxProblem), stations(station).y(idxProblem), 'k.', 'MarkerSize', 15)
    plot(stations(station).x(idxProblem-1), stations(station).y(idxProblem-1), 'g.', 'MarkerSize', 15)
    xlabel("X [km]"); ylabel("Y [km]")
    legend([sat, stat], ["Satellite", "Station " + num2str(station)])
    
    phi = stations(station).phi(idxProblem);
    x = X(1, idxProblem);
    y = X(3, idxProblem);
    y_s = stations(station).y(idxProblem);
    theta_s = stations(station).theta(idxProblem);

    % if (y-y_s) < 0 && y > 0 && x < 0 % Quadrant 2
    %     phi_check = phi + 2*pi;
    % elseif (y-y_s) > 0 && y < 0 && x < 0 % Quadrant 3
    %     phi_check = phi - 2*pi;
    % else
    %     phi_check = phi;
    % end
    % 
    % viewingAngle = phi_check - theta_s;
    viewingAngle = wrapTo2Pi(theta_s - phi);

    fprintf("station %.0f: visible: %.0f, view angle = %.3f, phi = %.3f, theta_s = %.3f\n", station, stations(station).visible(idxProblem), viewingAngle, stations(station).phi(idxProblem), stations(station).theta(idxProblem))

end

