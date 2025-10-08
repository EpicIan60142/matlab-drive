%% ASEN 6014 Priject 1 Main Script
% By: Ian Faber, 10/07/2025

%% Housekeeping
clc; clear; close all;

addpath("Plotting/")
addpath("RaceCourse/")
addpath("Dynamics/")
addpath("../utilities/")

%% Setup
    % Planetary constants structure
        % Point mass
pConst.mu = 398600.4415; % km^3/m^2
        % J2
pConst.J2 = 1.08264e-3;
pConst.Ri = 6378; % km
        % Drag
pConst.rho0 = 3.614e-22; % kg/km^3
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
T = 120*60; % 180 minute orbit period
oe_c.mu = pConst.mu;
oe_c.a = (pConst.mu*(T/(2*pi))^2)^(1/3);
oe_c.e = 0.1; % Decently eccentric course origin orbit
oe_c.i = deg2rad(45);
oe_c.RAAN = deg2rad(15);
oe_c.argPeri = deg2rad(20);
oe_c.f = deg2rad(25);

    % ode45 settings
opt = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);

%% Task 1: Generate random race course
    % Determine if a random or seeded course is made
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
courseFig = plotCourse(rings, figNum, titleText, xLabel, yLabel, zLabel, [], []);

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

%% Task 3/4: Develop integration scheme for all rings with J2 and drag and propagate over an average race
    % Compile initial ring states together
cart_c = convClassicOE2Cart(oe_c);
X_c = [cart_c.rVec; cart_c.vVec];
X0 = X_c;
for k = 1:length(rings)
    ringRho = [rings(k).center; zeros(3,1)]; % All rings start at their center positions at rest
    ringX0 = convDeputyH2N(X_c, ringRho, pConst);
    X0 = [X0; ringX0];
end

    % Turn on disturbances for each ring
J2_enable = true;
drag_enable = false;
disturb = [false; false]; % Course origin isn't subject to J2 or drag
for k = 1:length(rings)
    disturb = [disturb, [J2_enable; drag_enable]];
end

    % Define timespan
tspan = 0:1:17*60; % Average race of 17 minutes

    % Integrate ring motion forward in time in inertial coordinates
[t_Inertial, X_Inertial] = ode45(@(t,X)multiSatOrbitEOMInertial(t,X,pConst,[rings(1);rings],disturb), tspan, X0, opt);

    % Assign ring trajectories in the Hill frame
rings_Inertial = rings;
for k = 1:length(rings_Inertial)
    ringX = [];
    for kk = 1:length(t_Inertial)
        ringX = [ringX, convDeputyN2H(X_Inertial(kk,1:6)', X_Inertial(kk,(k*6 + 1):(k*6 + 6))',pConst)];
    end
    rings_Inertial(k).X = ringX;
    rings_Inertial(k).t = t_Inertial;
end

    % Plot numerically integrated rings
figNum = 3; 
titleText = sprintf("Generated Race Course with Ring Trajectories");
xLabel = sprintf("Radial [km]"); yLabel = sprintf("Along-Track [km]"); zLabel = sprintf("Cross-Track [km]");
trajStyle = 'b-'; trajLabel = sprintf("Ring Trajectories - Inertial mapped to Hill");
courseFig2 = plotCourse(rings_Inertial, figNum, titleText, xLabel, yLabel, zLabel, trajStyle, trajLabel);

    % Integrate ring motion forward in time in Hill frame
X0 = X_c;
for k = 1:length(rings)
    X0 = [X0; [rings(k).center; zeros(3,1)]];
end
[t_Hill, X_Hill] = ode45(@(t,X)multiSatOrbitEOMRelative(t,X,pConst,[rings(1);rings],disturb), tspan, X0, opt);

    % Assign ring trajectories in the Hill frame
rings_Hill = rings;
for k = 1:length(rings_Hill)
    rings_Hill(k).X = X_Hill(:, (6*k + 1):(6*k + 6))';
    rings_Hill(k).t = t_Hill;
end

    % Plot numerically integrated rings
figNum = 4; 
titleText = sprintf("Generated Race Course with Ring Trajectories");
xLabel = sprintf("Radial [km]"); yLabel = sprintf("Along-Track [km]"); zLabel = sprintf("Cross-Track [km]");
trajStyle = 'r-'; trajLabel = sprintf("Ring Trajectories - Integrated in Hill");
courseFig3 = plotCourse(rings_Hill, figNum, titleText, xLabel, yLabel, zLabel, trajStyle, trajLabel);

    % Check answer
figure; tl = tiledlayout(3,1);
nexttile;
    hold on; grid on;
    plot(t_Inertial, rings_Inertial(1).X(1,:) - rings_Hill(1).X(1,:))
nexttile;
    hold on; grid on;
    plot(t_Inertial, rings_Inertial(1).X(2,:) - rings_Hill(1).X(2,:))
nexttile;
    hold on; grid on;
    plot(t_Inertial, rings_Inertial(1).X(3,:) - rings_Hill(1).X(3,:))





