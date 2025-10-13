%% ASEN 6014 HW 2 Problem 2 Main Script
% By: Ian Faber, 10/03/2025

%% Housekeeping
clc; clear; close all;

%% Setup
    % Planetary constants
pConst.mu = 398600.4415; % km^3/s^2

    % Chief orbit
oe_c.mu = pConst.mu;
oe_c.a = 10000; % km
oe_c.e = 0.3; % <- Required by the problem!
oe_c.i = deg2rad(15);
oe_c.RAAN = deg2rad(69);
oe_c.argPeri = deg2rad(42);
oe_c.f = deg2rad(5);

%         % Logan params
% oe_c.a = 8000;
% oe_c.i = deg2rad(5);
% oe_c.RAAN = deg2rad(0);
% oe_c.argPeri = deg2rad(0);
% oe_c.f = deg2rad(0);

    % Bounded formation requirements
zBound = 0.25; % km, max out of plane oscillation
xyBound = 1; % km, max in plane oscillation
yGoal = 0.5; % km, desired along track offset

    % Intermediate parameters
eta = sqrt(1-oe_c.e^2);

%% Calculate orbit element differences for a specific di and dargPeri for 
% the given desired bounded formation
di = linspace(1e-3, -1e-3, 10);
dargPeri = linspace(1e-3, -1e-3, 10);

for k = 1:length(di)
        % Choose di and dargPeri
    doe.di = deg2rad(di(k));
    doe.dargPeri = deg2rad(dargPeri(k));
    
        % Calculate other elements based on di and dargPeri
    doe.da = 0; % Required for bounded motion
    doe.dRAAN = sqrt(((zBound/(oe_c.a*(1+oe_c.e)))^2 - doe.di^2)/(sin(oe_c.i)^2));
    doe.dM0 = ((2*eta^3)/(2+oe_c.e))*((yGoal/(oe_c.a*(1+oe_c.e))) - doe.dargPeri - cos(oe_c.i)*doe.dRAAN);
    doe.de = sqrt(((2*(eta^2)*xyBound)/(oe_c.a*(1+oe_c.e)*(4+oe_c.e)))^2 - ((oe_c.e^2)*doe.dM0^2)/(eta^2));
    
        % Populate deputy spacecraft orbit elements
    oe_d.mu = pConst.mu;
    oe_d.a = oe_c.a + doe.da;
    oe_d.e = oe_c.e + doe.de;
    oe_d.i = oe_c.i + doe.di;
    oe_d.RAAN = oe_c.RAAN + doe.dRAAN;
    oe_d.argPeri = oe_c.argPeri + doe.dargPeri;
    M_c = convE2M(convf2E(oe_c.f, oe_c.e), oe_c.e);
    M_d = M_c + doe.dM0;
    E_d = convM2E(M_d, oe_d.e, false);
    oe_d.f = convE2f(E_d, oe_d.e);
    
    %% Simulate via linearized analytical solution
        % True anomalies for simulation
    f = linspace(0, 2*pi, 100);
    
        % Calculate analytical solution
    rho_lin = calcLinearGenEllSol(oe_c, doe, f);
    
    %% Nuimerically integrate
        % Find orbit period
    T = 2*pi*sqrt((oe_c.a^3)/pConst.mu);
    tspan = 0:1:2*T;
    
        % Collate start states
    cart_c = convClassicOE2Cart(oe_c);
    X_c = [cart_c.rVec; cart_c.vVec];
    
    cart_d = convClassicOE2Cart(oe_d);
    X_d = [cart_d.rVec; cart_d.vVec];
    X_d = convDeputyN2H(X_c, X_d, pConst);
    
    X0 = [X_c; X_d];
    
        % Set ODE options
    opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
    
        % Integrate
    [t_NL, X_NL] = ode45(@(t,X)fullNLGenRelEOM(t,X,pConst), tspan, X0, opt);
    
    %% Plot solutions
        % Plot 3D solutions
    fig = figure;
    hold on; grid on;
    title("Bounded Formation")
    plot3(0,0,0,'k*')
    plot3(0,yGoal,0,'m*')
    plot3(rho_lin(1,:), rho_lin(2,:), rho_lin(3,:), 'b-');
    plot3(X_NL(:,7), X_NL(:,8), X_NL(:,9), 'r--');
    xlabel("Radial [km]"); ylabel("Along-Track [km]"); zlabel("Cross-Track [km]")
    legend("Chief spacecraft", "Desired Offset", "Linearized", "Full NL")
    view([30 35])
    
    if k == 1
        %% Plot vs. time
    figure; tl = tiledlayout(3,1); ax = [];
    title(tl, "Nonlinear Deputy Position vs. Time for 2 Orbit Periods")
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        title("x vs. time")
        plot(t_NL, X_NL(:,7), 'b-');
        yline(xyBound, 'k--');
        yline(-xyBound, 'k--');
        xlabel("Time [sec]"); ylabel("x [km]");
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        title("y vs. time")
        plot(t_NL, X_NL(:,8), 'b-');
        yline(yGoal+xyBound, 'k--');
        yline(yGoal-xyBound, 'k--');
        a = yline(yGoal, 'm--');
        xlabel("Time [sec]"); ylabel("y [km]");
        legend(a, "Desired Along-Track Offset", 'Location', 'best');
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        title("z vs. time")
        plot(t_NL, X_NL(:,9), 'b-');
        yline(zBound, 'k--');
        yline(-zBound, 'k--');
        xlabel("Time [sec]"); ylabel("z [km]");
    end
end

