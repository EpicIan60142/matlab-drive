%% ASEN 6014 HW 3 Problem 3 Main Script
% By: Ian Faber, 11/08/2025

%% Housekeeping
clc; clear; close all;

addpath("../utilities/")

%% Setup
    % Planetary constants structure
pConst.mu = 398600.4415; % km^3/s^2
pConst.Ri = 6378; % km
pConst.J2 = 1.08264e-3;

    % Chief initial orbit elements
oe_c0.mu = pConst.mu;
oe_c0.a = 8500; % km
oe_c0.e = 0.05;
oe_c0.i = deg2rad(53); % deg -> rad
oe_c0.RAAN = deg2rad(55); % deg -> rad
oe_c0.argPeri = deg2rad(40); % deg -> rad
oe_c0.f = deg2rad(0); % deg -> rad

    % Initial deputy orbit element differences
doe.da = 1; % km
doe.de = 1e-3;
doe.di = deg2rad(0.1);
doe.dRAAN = deg2rad(0.1);

    % Initial deputy orbit elements
oe_d0 = oe_c0;
oe_d0.a = oe_c0.a + doe.da;
oe_d0.e = oe_c0.e + doe.de;
oe_d0.i = oe_c0.i + doe.di;
oe_d0.RAAN = oe_c0.RAAN + doe.dRAAN;

    % Desired orbit element differences
doe_r.da = 0*5e-2;
doe_r.de = 1*5e-4;
doe_r.di = 0*deg2rad(0.1);
doe_r.dRAAN = 0*deg2rad(0.1);
doe_r.dargPeri = 0*deg2rad(0.1);
doe_r.dM0 = 0*deg2rad(0.01);

    % Controller parameters
kConst.K = 1e-10;
kConst.P11 = 0.004;
kConst.P11Off = 1e-12;
kConst.P11Amp = 1e-3;
kConst.P11Exp = 10;
kConst.P22 = 0.008;
kConst.P22Off = 1e-12;
kConst.P22Amp = 1e-3;
kConst.P22Exp = 10;
kConst.P33 = 0.004;
kConst.P33Off = 1e-12;
kConst.P33Amp = 1e-3;
kConst.P33Exp = 10;
kConst.P44Off = 1e-12;
kConst.P44Amp = 1e-3;
kConst.P44Exp = 10;
kConst.P55Off = 1e-12;
kConst.P55Amp = 1e-3;
kConst.P55Exp = 10;
kConst.P66Off = 1e-12;
kConst.P66Amp = 1e-3;
kConst.P66Exp = 10;

    % ode45 settings
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);

%% Run controller
    % Setup initial state vector
cart_c = convClassicOE2Cart(oe_c0);
X0_c = [cart_c.rVec; cart_c.vVec];

cart_d = convClassicOE2Cart(oe_d0);
X0_d = [cart_d.rVec; cart_d.vVec];

X0 = [X0_c; X0_d];

    % Define tspan for multiple orbits
nOrbits = 20;
T = 2*pi*sqrt((oe_c0.a)^3/pConst.mu);
tspan = 0:10:nOrbits*T;

    % Run controller
J2Dist = true;
% controller = @(t,X) OEFeedbackControlElems(t,X,doe_r,kConst,pConst,J2Dist);
controller = @(t,X) OEAltFeedbackControlElems(t,X,doe_r,kConst,pConst,J2Dist);
[t, X] = ode45(@(t,X)controller(t,X), tspan, X0, opt);

[~, u, doe, oe_d, oe_c] = cellfun(@(t,X)controller(t,X.'), num2cell(t), num2cell(X,2), 'uni', 0);
u = cellfun(@(x)x',u,'uni',0);
doe = cellfun(@(x)x',doe,'uni',0);
oe_d = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_d, 'uni', 0);
oe_c = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_c, 'uni', 0);

%% Plot and report results
    % Extract states and control
X_c = X(:,1:6);
X_d = X(:,7:12);
u = cell2mat(u);
doe = cell2mat(doe);
oe_d = cell2mat(oe_d);
oe_c = cell2mat(oe_c);

    % Convert to Hill Frame
X_d_Hill = zeros(size(X_d));
for k = 1:length(t)
    X_d_Hill(k,:) = convDeputyN2H(X_c(k,:)', X_d(k,:)', pConst)';
end

    % Plot
markerSize = 30;
figure; axis equal
hold on; grid on;
title("OE Difference Alternate Relative Motion Control")
scatter3(X_d_Hill(1,1), X_d_Hill(1,2), X_d_Hill(1,3), markerSize, 'g', 'filled');
plot3(X_d_Hill(:,1), X_d_Hill(:,2), X_d_Hill(:,3), 'b-')
scatter3(X_d_Hill(end,1), X_d_Hill(end,2), X_d_Hill(end,3), markerSize, 'r', 'filled')
% scatter3(0,0,0,markerSize,'k','*')
xlabel("Radial [km]"); ylabel("Along-Track [km]"); zlabel("Cross-Track [km]")
view([30,35]);
legend("Start", "Trajectory", "End")%, "Chief")

figure; tl = tiledlayout(3,1); ax = [];
title(tl, "Control Effort vs. Time")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t,u(:,1),'b-');
    xlabel("Time [sec]"); ylabel("u_1 [m/s^2]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t,u(:,2),'b-');
    xlabel("Time [sec]"); ylabel("u_2 [m/s^2]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t,u(:,3),'b-');
    xlabel("Time [sec]"); ylabel("u_3 [m/s^2]");
linkaxes(ax, 'x');

figure; tl = tiledlayout('flow'); ax = [];
title(tl, "Difference in Actual and Desired Deputy Elements vs. Time")
try doe_r.da;
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        plot(t,doe(:,1),'b-');
        plot(t,zeros(size(t)), 'r--')
        xlabel("Time [sec]"); ylabel("\Delta a/r_e")
end
try doe_r.de;
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        plot(t,doe(:,2),'b-');
        plot(t,zeros(size(t)), 'r--')
        xlabel("Time [sec]"); ylabel("\Delta e")
end
try doe_r.di;
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        plot(t,doe(:,3),'b-');
        plot(t,zeros(size(t)), 'r--')
        xlabel("Time [sec]"); ylabel("\Delta i [rad]")
end
try doe_r.dRAAN;
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        plot(t,doe(:,4),'b-');
        plot(t,zeros(size(t)), 'r--')
        xlabel("Time [sec]"); ylabel("\Delta \Omega [rad]")
end
try doe_r.dargPeri;
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        plot(t,doe(:,5),'b-');
        plot(t,zeros(size(t)), 'r--')
        xlabel("Time [sec]"); ylabel("\Delta \omega [rad]")
end
try doe_r.dM0;
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        plot(t,doe(:,6),'b-');
        plot(t,zeros(size(t)), 'r--')
        xlabel("Time [sec]"); ylabel("\Delta M [rad]")
end
linkaxes(ax, 'x');

figure; tl = tiledlayout(3,2); ax = [];
title(tl, "Difference in Deputy and Chief Elements vs. time")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t,oe_d(:,1) - oe_c(:,1),'b-');
    try doe_r.da; 
        plot(t,doe_r.da*ones(size(t)), 'r--') 
    end
    xlabel("Time [sec]"); ylabel("a_d - a_c [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t,oe_d(:,2) - oe_c(:,2),'b-');
    try doe_r.de;
        plot(t,doe_r.de*ones(size(t)), 'r--')
    end
    xlabel("Time [sec]"); ylabel("e_d - e_c")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t,oe_d(:,3) - oe_c(:,3),'b-');
    try doe_r.di;
        plot(t,doe_r.di*ones(size(t)), 'r--')
    end
    xlabel("Time [sec]"); ylabel("i_d - i_c [rad]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t,oe_d(:,4) - oe_c(:,4),'b-');
    try doe_r.dRAAN;
        plot(t,doe_r.dRAAN*ones(size(t)), 'r--')
    end
    xlabel("Time [sec]"); ylabel("\Omega_d - \Omega_c [rad]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t,oe_d(:,5) - oe_c(:,5),'b-');
    try doe_r.dargPeri;
        plot(t,doe_r.dargPeri*ones(size(t)), 'r--')
    end
    xlabel("Time [sec]"); ylabel("\omega_d - \omega_c [rad]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t,oe_d(:,6) - oe_c(:,6),'b-');
    try doe_r.dM0;
        plot(t,doe_r.dM0*ones(size(t)), 'r--')
    end
    xlabel("Time [sec]"); ylabel("M_d - M_c [rad]")
linkaxes(ax, 'x');

frames = 10;
figure;
for k = 1:frames:length(t)
    clf;
    hold on; grid on; axis equal
    title(sprintf("Deputy Hill Trajectory at t = %.3f sec", t(k)));
    cutoff = find(t >= 3*T, 1, 'first');
    if k <= cutoff
        start = 1;
        plot3(X_d_Hill(1,1), X_d_Hill(1,2), X_d_Hill(1,3), 'g.', 'MarkerSize', markerSize);
    else
        start = floor(k - cutoff);
    end
    plot3(X_d_Hill(start:k,1), X_d_Hill(start:k,2), X_d_Hill(start:k,3), 'b-');
    if k <= length(t) - frames
        plot3(X_d_Hill(k,1), X_d_Hill(k,2), X_d_Hill(k,3), 'k.', 'MarkerSize', markerSize);
    end
    plot3(0,0,0,'k*','MarkerSize',10);
    if k >= length(t) - frames
       plot3(X_d_Hill(end,1), X_d_Hill(end,2), X_d_Hill(end,3), 'r.', 'MarkerSize', markerSize); 
    end
    xlabel("Radial [km]"); ylabel("Along-Track [km]"); zlabel("Cross-Track [km]");
    view([30,35]);
    drawnow;
end



