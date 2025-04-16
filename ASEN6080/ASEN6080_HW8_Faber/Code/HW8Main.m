%% ASEN 6080 HW 8 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
    % Path logistics
addpath("..\")
addpath(genpath("..\..\..\Utilities\"));

    % Extract Earth params
pConst = getPlanetConst();

    % Extract spacecraft params
scConst = getSCConst();

    % Extract orbit params
orbital = getOrbitConst();

    % Number of Monte Carlo runs and time to propagate each run for
N = 3000;
filename = sprintf("..\\Data\\MonteCarlo%.0fRuns.mat", N);
tspan = 0:10:24*60*60; % 24 hours in sec
titleText = sprintf("%.0f Monte Carlo Orbits", N);

    % Covariance
sigR = 1; % km
sigV = 0.001; % km/s
P0 = diag([sigR^2, sigR^2, sigR^2, sigV^2, sigV^2, sigV^2]);

    % Earth plotting shenanigans
theta0 = 0;
[earthX_s, earthY_s, earthZ_s] = sphere(40);
earthX_s = pConst.Ri*earthX_s;
earthY_s = pConst.Ri*earthY_s;
earthX = earthX_s*cos(theta0) - earthY_s*sin(theta0);
earthY = earthX_s*sin(theta0) + earthY_s*cos(theta0);
earthZ = pConst.Ri*earthZ_s;
earthMap = imread('2k_earth_daymap.jpg'); % Image from https://www.solarsystemscope.com/textures/

%% Part 1a/b. Generate Monte Carlo orbit data
rng(69420); % Set rng seed for consistency

    % Generate new data by changing the boolean below to true
if false
    runs = [];
    for k = 1:N
        x0Perturb = chol(P0)*randn(6,1);
    
        data = generateTruthOrbit_MuJ2Drag(pConst, orbital, scConst, x0Perturb, [], tspan);
    
        runs = [runs; data];
    
        fprintf("Generated %.0f orbits so far\n", k)
    end
    
    save(filename, "runs", '-mat');
else
    load("..\Data\MonteCarlo3000Runs.mat");
    fprintf("Loaded Data\n")
end

    % Create nominal orbit
x0Perturb = zeros(6,1);
nominal = generateTruthOrbit_MuJ2Drag(pConst, orbital, scConst, x0Perturb, [], tspan);

    % Plot orbits
figure
hold on; grid on; axis equal;
title(titleText);
earth = surf(earthX, earthY, earthZ, 'FaceAlpha', 0.5);
set(earth,'FaceColor','texturemap','cdata',earthMap,'edgecolor','none');
for k = 1:length(runs)
    plot3(runs(k).X_ref(:,1), runs(k).X_ref(:,2), runs(k).X_ref(:,3))
end
scale = 1.25;
xlim([-scale*orbital.a, scale*orbital.a]); ylim([-scale*orbital.a, scale*orbital.a]); zlim([-scale*orbital.a, scale*orbital.a]);
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
view([30 35]);

%% Part 1c. Analyze runs at 6 hour intervals
analysisTimes = linspace(0,24,5)*60*60; % Convert to sec

nomTraj = [];
meanTraj = [];
stdTraj = [];
PTraj = {};
for k = 1:length(analysisTimes)
        % Find the index corresponding to the desired analysis time
    idx = find(tspan == analysisTimes(k));

        % Pull out nominal trajectory at this time
    monteTraj = [];
    for kk = 1:length(runs)
        monteTraj = [monteTraj, runs(kk).X_ref(idx,:)'];
    end
    nom = nominal.X_ref(idx,:)';

        % Calculate mean and standard deviation of trajectories
    mu = mean(monteTraj,2);
    sigma = std(monteTraj,0,2);

        % Construct covariance matrix
    P = cov(monteTraj'); % Need input to be structured with observations (runs) as rows and variables (X, Y, etc.) as columns

        % Make corner plot
    cornerPlot(nom, monteTraj, analysisTimes(k));

        % Save nominal, mean, std, and covariance matrix
    nomTraj = [nomTraj, nom];
    meanTraj = [meanTraj, mu];
    stdTraj = [stdTraj, sigma];
    PTraj = [PTraj, {P}];

end





