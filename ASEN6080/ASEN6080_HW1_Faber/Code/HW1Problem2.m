%% ASEN 6080 HW 1 Problem 2 Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Part 2a
    % Gravitational parameters
mu = 398600.4415; % km^3/s^2
J2 = 1.08264e-3; 
Ri = 6378; % km

    % Orbit parameters
a = 10000; % km
e = 0.001; 
i = deg2rad(40); % deg -> rad
RAAN = deg2rad(80); % deg -> rad
argPeri = deg2rad(40); % deg -> rad
truAnom = deg2rad(0); % deg -> rad

    % Convert to cartesian and build initial state
orbital = [mu; a; e; i; RAAN; argPeri; truAnom];
X0_nom = convOrbitalToCartesian(orbital);
X0_nom = [X0_nom; J2];

    % Propagate for 15 orbits
T = 2*pi*sqrt(a^3/mu);
dt = T/1000;

tspan = 0:dt:(15*T);
opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

[t, X_nom] = ode45(@(t,X)orbitEOM_MuJ2(t,X,mu,Ri), tspan, X0_nom, opt);

    % Plot trajectory
figure
hold on; grid on; axis equal;
title("Orbit Trajectory - Nominal");
plot3(X_nom(:,1), X_nom(:,2), X_nom(:,3));
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]")
view([30 35]);

    % Add small perturbation to initial state and propagate again
dx0 = [1; 0; 0; 0; 10e-3; 0; 0];
X0_perturb = X0_nom + dx0;

[t_perturb, X_perturb] = ode45(@(t,X)orbitEOM_MuJ2(t,X,mu,Ri), tspan, X0_perturb, opt);

    % Plot trajectory
figure
hold on; grid on; axis equal;
title("Orbit Trajectory - Perturbed");
plot3(X_perturb(:,1), X_perturb(:,2), X_perturb(:,3));
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]")
view([30 35]);

    % Extract perturbations
dx = X_perturb - X_nom;

%% Part 2b
X0(1) = -0.85188696962247;
X0(2) = 0.80032070980182;
X0(3) = -1.50940472473439;
X0(4) = 0.87587414783453;
X0(5) = -0.24278953633334;
X0(6) = 0.16681343945350;
X0(7) = -1.96541870928278;
Phi0(1,1) = -1.27007139263854;
Phi0(1,2) = -1.86512257453063;
Phi0(1,3) = 0.06600934128821;
Phi0(1,4) = 0.59069655120545;
Phi0(1,5) = -0.34563190830705;
Phi0(1,6) = -1.47629235201010;
Phi0(1,7) = 0.59430761682985;
Phi0(2,1) = 1.17517126546302;
Phi0(2,2) = -1.05110705924059;
Phi0(2,3) = 0.45129021363078;
Phi0(2,4) = -0.63578573784723;
Phi0(2,5) = -1.17140482049761;
Phi0(2,6) = 0.25889995716040;
Phi0(2,7) = -0.27646490663926;
Phi0(3,1) = 2.02916018474976;
Phi0(3,2) = -0.41738204799680;
Phi0(3,3) = -0.32220971801190;
Phi0(3,4) = 0.60334661284576;
Phi0(3,5) = -0.68558678043728;
Phi0(3,6) = -2.01869095243834;
Phi0(3,7) = -1.85758288592737;
Phi0(4,1) = -0.27515724067569;
Phi0(4,2) = 1.40216228633781;
Phi0(4,3) = 0.78840921622743;
Phi0(4,4) = -0.53524796777590;
Phi0(4,5) = 0.92621639416896;
Phi0(4,6) = 0.19974026229838;
Phi0(4,7) = 0.04073081174943;
Phi0(5,1) = 0.60365844582581;
Phi0(5,2) = -1.36774699097611;
Phi0(5,3) = 0.92873604681331;
Phi0(5,4) = -0.15508038549279;
Phi0(5,5) = -1.48167521167231;
Phi0(5,6) = 0.42586431913121;
Phi0(5,7) = 0.28297017716199;
Phi0(6,1) = 1.78125189324250;
Phi0(6,2) = -0.29253499915187;
Phi0(6,3) = -0.49079037626976;
Phi0(6,4) = 0.61212237077216;
Phi0(6,5) = -0.55805780868504;
Phi0(6,6) = -1.27004345059705;
Phi0(6,7) = 0.06356121930250;
Phi0(7,1) = 1.77365832632615;
Phi0(7,2) = 1.27084843418894;
Phi0(7,3) = 1.79720058425494;
Phi0(7,4) = -1.04434349451734;
Phi0(7,5) = -0.02845311157066;
Phi0(7,6) = -0.48521883574304;
Phi0(7,7) = 0.43343006511160;

XPhi0 = [reshape(X0,7,1); reshape(Phi0,49,1)];

dXPhi = STMEOM_J2(0,XPhi0,mu,Ri);

XDot = dXPhi(1:7)
PhiDot_mat = reshape(dXPhi(8:56),7,7)

%% Part 2c
    % Integrate Phi across entire timespan
Phi0_nom = eye(7);

XPhi0 = [X0_nom; reshape(Phi0_nom,49,1)];

[t_phi, XPhi] = ode45(@(t,XPhi)STMEOM_J2(t,XPhi,mu,Ri), tspan, XPhi0, opt);

    % Propagate perturbations
dx_phi = [];
for k = 1:size(XPhi,1)
    Phi = reshape(XPhi(k,8:56),7,7);
    dx_phi = [dx_phi; (Phi*dx0)'];
end

    % Plot perturbations
figure; t = tiledlayout(4,2);
title(t, "\deltaX vs. time")
nexttile
    hold on; grid on;
    plot(t_perturb, dx(:,1),'b-')
    plot(t_phi, dx_phi(:,1),'r--')
    xlabel("Time [sec]"); ylabel("\deltax [km]");
nexttile
    hold on; grid on;
    method1 = plot(t_perturb, dx(:,4),'b-');
    method2 = plot(t_phi, dx_phi(:,4),'r--');
    xlabel("Time [sec]"); ylabel("\deltaxDot [km/s]");
nexttile
    hold on; grid on;
    plot(t_perturb, dx(:,2),'b-')
    plot(t_phi, dx_phi(:,2),'r--')
    xlabel("Time [sec]"); ylabel("\deltay [km]");
nexttile
    hold on; grid on;
    plot(t_perturb, dx(:,5),'b-')
    plot(t_phi, dx_phi(:,5),'r--')
    xlabel("Time [sec]"); ylabel("\deltayDot [km/s]");
nexttile
    hold on; grid on;
    plot(t_perturb, dx(:,3),'b-')
    plot(t_phi, dx_phi(:,3),'r--')
    xlabel("Time [sec]"); ylabel("\deltaz [km]");
nexttile
    hold on; grid on;
    plot(t_perturb, dx(:,6),'b-')
    plot(t_phi, dx_phi(:,6),'r--')
    xlabel("Time [sec]"); ylabel("\deltazDot [km/s]");
nexttile
    hold on; grid on;
    plot(t_perturb, dx(:,7),'b-');
    plot(t_phi, dx_phi(:,7),'r--');
    xlabel("Time [sec]"); ylabel("\deltaJ_2");
legend([method1, method2], ["Integrator Method", "STM Method"], 'Location', 'bestoutside');

%% Part 2d
    % Find difference vector
dx_diff = dx - dx_phi;

    % Plot difference vector
figure; t = tiledlayout(4,2);
title(t, "(\deltaX_{Int} - \deltaX_{\Phi}) vs. time")
nexttile
    hold on; grid on;
    plot(t_phi, dx_diff(:,1),'b-')
    xlabel("Time [sec]"); ylabel("\deltax Diff [km]")
nexttile
    hold on; grid on;
    plot(t_phi, dx_diff(:,4),'b-')
    xlabel("Time [sec]"); ylabel("\deltaxDot Diff [km/s]")
nexttile
    hold on; grid on;
    plot(t_phi, dx_diff(:,2),'b-')
    xlabel("Time [sec]"); ylabel("\deltay Diff [km]")
nexttile
    hold on; grid on;
    plot(t_phi, dx_diff(:,5),'b-')
    xlabel("Time [sec]"); ylabel("\deltayDot Diff [km/s]")
nexttile
    hold on; grid on;
    plot(t_phi, dx_diff(:,3),'b-')
    xlabel("Time [sec]"); ylabel("\deltaz Diff [km]")
nexttile
    hold on; grid on;
    plot(t_phi, dx_diff(:,6),'b-')
    xlabel("Time [sec]"); ylabel("\deltazDot Diff [km/s]")
nexttile
    hold on; grid on;
    plot(t_phi, dx_diff(:,7),'b-')
    xlabel("Time [sec]"); ylabel("\deltaJ_2 Diff")
