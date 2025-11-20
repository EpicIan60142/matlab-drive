%% ASEN 6014 Final Project Main Script
% By: Ian Faber, 11/19/2025

%% Housekeeping
clc; clear; close all;

addpath("Controls\")
addpath("Dynamics\")
addpath("Plotting\")
addpath("RaceCourse\")
addpath("../utilities/")

%% Setup
    % Planetary constants structure
        % Point mass
pConst.mu = 398600.4415; % km^3/m^2
        % J2
pConst.J2 = 1.08264e-3;
pConst.Ri = 6378; % km
        % Drag
pConst.rho0 = 3.614e-4; % kg/km^3
pConst.r0 = 700 + pConst.Ri; % km
pConst.H = 88.667; % km

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

        % Recover number of rings
numRings = numRings + 2;

        % Course origin orbit settings
T = 180*60; % 180 minute orbit period
oe_c.mu = pConst.mu;
oe_c.a = (pConst.mu*(T/(2*pi))^2)^(1/3);
oe_c.e = 0.05; % Decently eccentric course origin orbit
oe_c.i = deg2rad(45);
oe_c.RAAN = deg2rad(15);
oe_c.argPeri = deg2rad(-23);
oe_c.f = deg2rad(10);

    % ode45 settings
opt = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);

    % Video saving
videoFolder = "Videos/";
%% Task 0: Generate a desired race course formation
    % Determine if a random course is made
randomCourse = false;

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
startRing.m = 100;
rings(1).params.lastRing = startRing;
rings = [startRing; rings];

        % End ring
endRing = generateRing(min(semiMaj), min(semiMin), 0, 0, min(interRingDist), rings(end));
endRing.m = 1;
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

%% Task 1a: Design and implement a race course formation feedback control law
    % Define gains, goal is to go from lead-follower to end position within
    % 10 minutes
kConst_course.K1 = 2e-5;
kConst_course.K2 = 5e-3;

    % Get initial state of the race course origin
cart_c = convClassicOE2Cart(oe_c);
X_c0 = [cart_c.rVec; cart_c.vVec];

    % Test ring formation control on a random ring from the race course
X_d0 = convDeputyH2N(X_c0, zeros(6,1), pConst); % Start at the course origin
X0_test = [X_c0; X_d0];
tspan_test = 0:10:T; % Ideally want rings to reach their positions within an hour
ringIdx = randi(length(rings),1);

[t_test, X_test] = ode23t(@(t,X)ringFormationFeedbackControl(t, X, rings(ringIdx).center, kConst_course, pConst, rings(ringIdx)), tspan_test, X0_test);
[~, u_test] = cellfun(@(t,X)ringFormationFeedbackControl(t, X.', rings(ringIdx).center, kConst_course, pConst, rings(ringIdx)), num2cell(t_test), num2cell(X_test,2), 'uni', 0);
u_test = cellfun(@(x)x', u_test, 'UniformOutput', false);
u_test = cell2mat(u_test);

    % Convert deputy from Inertial frame to Hill frame
X_test_Hill = [];
for k = 1:length(t_test)
    X_test_Hill = [X_test_Hill; convDeputyN2H(X_test(k, 1:6)', X_test(k, 7:12)', pConst)'];
end

    % Assign trajectory arrays
rings(ringIdx).X = X_test_Hill;
rings(ringIdx).t = t_test;
rings(ringIdx).u = u_test;
rings(ringIdx).X_r = repmat([rings(ringIdx).center; zeros(3,1)]', length(t_test), 1);
rings(ringIdx).oe = NaN(length(t_test), 6); % Not controlling orbit elements
rings(ringIdx).oe_r = NaN(length(t_test), 6);

%% Task 1b: Design and implement a lead-follower formation feedback control law
    % Define lead follower characteristics
doe_r.da = 0; % km
doe_r.de = 0;
doe_r.di = deg2rad(0); % deg -> rad
doe_r.dRAAN = deg2rad(0); % deg -> rad
doe_r.dargPeri = deg2rad(0); % deg -> rad
doe_r.dM0 = deg2rad(2e-4); % deg -> rad

kConst_lead.P11 = 4e-3;
kConst_lead.P22 = 8e-3;
kConst_lead.P33 = 4e-3;

    % Add another orbit of time
tspan_test = (0:10:3*T) + t_test(end); % Want rings to stow within two orbits

    % Update initial state
X0_test = X_test(end,:)';

    % Run lead-follower control
[t_test, X_test] = ode45(@(t,X)leadFollowerFeedbackControl(t, X, doe_r, kConst_lead, pConst, rings(ringIdx)), tspan_test, X0_test, opt);
[~, u_test, doe_test, oe_d_test, oe_c_test, oe_r_test] = cellfun(@(t,X)leadFollowerFeedbackControl(t, X.', doe_r, kConst_lead, pConst, rings(ringIdx)), num2cell(t_test), num2cell(X_test,2), 'uni', 0);
u_test = cellfun(@(x)x', u_test, 'UniformOutput', false);
u_test = cell2mat(u_test);
doe_test = cellfun(@(x)x',doe_test,'uni',0);
doe_test = cell2mat(doe_test);
oe_d_test = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_d_test, 'uni', 0);
oe_d_test = cell2mat(oe_d_test);
oe_c_test = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_c_test, 'uni', 0);
oe_c_test = cell2mat(oe_c_test);
oe_r_test = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_r_test, 'uni', 0);
oe_r_test = cell2mat(oe_r_test);

    % Convert deputy from Inertial frame to Hill frame
X_test_Hill = [];
for k = 1:length(t_test)
    X_test_Hill = [X_test_Hill; convDeputyN2H(X_test(k, 1:6)', X_test(k, 7:12)', pConst)'];
end

    % Assign trajectory arrays
rings(ringIdx).X = [rings(ringIdx).X; X_test_Hill];
rings(ringIdx).t = [rings(ringIdx).t; t_test];
rings(ringIdx).u = [rings(ringIdx).u; u_test];
rings(ringIdx).X_r = [rings(ringIdx).X_r; NaN(length(t_test), 6)]; % Not controlling Hill coords
rings(ringIdx).oe = [rings(ringIdx).oe; oe_d_test];
rings(ringIdx).oe_r = [rings(ringIdx).oe_r; oe_r_test];

%% Task 1c: Plot results of both controllers
    % Plot results
figNum = 1;
titleText = sprintf("Test of ring formation control");
xLabel = sprintf("Radial [km]"); yLabel = sprintf("Along-Track [km]"); zLabel = sprintf("Cross-Track [km]");
trajStyle = "b-"; trajLabel = sprintf("Ring Trajectory");
plotCourse(rings, figNum, titleText, xLabel, yLabel, zLabel, trajStyle, trajLabel);

    % Plot control effort
figure; tl = tiledlayout(3,1); ax = [];
title(tl, "Test Control Effort vs. Time");
nt = nexttile; ax = [ax; nt];
    semilogy(rings(ringIdx).t, abs(rings(ringIdx).u(:,1)), 'b-');
    hold on; grid on;
    xlabel("Time [sec]"); ylabel("u_r [km/s^2]");
nt = nexttile; ax = [ax; nt];  
    semilogy(rings(ringIdx).t, abs(rings(ringIdx).u(:,2)), 'b-');
    hold on; grid on;
    xlabel("Time [sec]"); ylabel("u_\theta [km/s^2]");
nt = nexttile; ax = [ax; nt];
    semilogy(rings(ringIdx).t, abs(rings(ringIdx).u(:,3)), 'b-');
    hold on; grid on;
    xlabel("Time [sec]"); ylabel("u_h [km/s^2]");
linkaxes(ax, 'x');

    % Plot ring trajectory
figure; tl = tiledlayout(3,2);
title(tl, "Test Ring Trajectory vs. Time")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).X(:,1), 'b-');
    plot(rings(ringIdx).t, rings(ringIdx).X_r(:,1), 'r--');
    xlabel("Time [sec]"); ylabel("x [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).X(:,4), 'b-');
    plot(rings(ringIdx).t, rings(ringIdx).X_r(:,4), 'r--');
    xlabel("Time [sec]"); ylabel("xDot [km/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).X(:,2), 'b-');
    plot(rings(ringIdx).t, rings(ringIdx).X_r(:,2), 'r--');
    xlabel("Time [sec]"); ylabel("y [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).X(:,5), 'b-');
    plot(rings(ringIdx).t, rings(ringIdx).X_r(:,5), 'r--');
    xlabel("Time [sec]"); ylabel("yDot [km/s]");
    legend("Trajectory", "Reference", 'location', 'eastoutside')
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).X(:,3), 'b-');
    plot(rings(ringIdx).t, rings(ringIdx).X_r(:,3), 'r--');
    xlabel("Time [sec]"); ylabel("z [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).X(:,6), 'b-');
    plot(rings(ringIdx).t, rings(ringIdx).X_r(:,6), 'r--');
    xlabel("Time [sec]"); ylabel("zDot [km/s]");
linkaxes(ax, 'x');

    % Plot ring elements
figure; tl = tiledlayout(3,2);
title(tl, "Test Ring Orbit Element Differences vs. Time")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).oe(:,1) - rings(ringIdx).oe_r(:,1), 'b-');
    plot(rings(ringIdx).t, 0*rings(ringIdx).oe_r(:,1), 'r--');
    xlabel("Time [sec]"); ylabel("\Delta a_d [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).oe(:,2) - rings(ringIdx).oe_r(:,2), 'b-');
    plot(rings(ringIdx).t, 0*rings(ringIdx).oe_r(:,2), 'r--');
    xlabel("Time [sec]"); ylabel("\Delta e_d");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rad2deg(rings(ringIdx).oe(:,3) - rings(ringIdx).oe_r(:,3)), 'b-');
    plot(rings(ringIdx).t, 0*rings(ringIdx).oe_r(:,3), 'r--');
    xlabel("Time [sec]"); ylabel("\Delta i_d [deg]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rad2deg(rings(ringIdx).oe(:,4) - rings(ringIdx).oe_r(:,4)), 'b-');
    plot(rings(ringIdx).t, 0*rings(ringIdx).oe_r(:,4), 'r--');
    xlabel("Time [sec]"); ylabel("\Delta \Omega_d [deg]");
    legend("Trajectory", "Reference", 'location', 'eastoutside')
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rad2deg(rings(ringIdx).oe(:,5) - rings(ringIdx).oe_r(:,5)), 'b-');
    plot(rings(ringIdx).t, 0*rings(ringIdx).oe_r(:,5), 'r--');
    xlabel("Time [sec]"); ylabel("\Delta \omega_d [deg]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rad2deg(rings(ringIdx).oe(:,6) - rings(ringIdx).oe_r(:,6)), 'b-');
    plot(rings(ringIdx).t, 0*rings(ringIdx).oe_r(:,6), 'r--');
    xlabel("Time [sec]"); ylabel("\Delta M_d [deg]");
linkaxes(ax, 'x');

    % Animate ring
animateRings(rings, true, videoFolder + "ControllerTest.mp4", 5, "");
