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
analysisTimes = linspace(0,24,5)*60*60; % Convert to sec

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

    % Dynamics function handle for UKF
DynFunc = @(t,X)orbitEOM_MuJ2Drag(t,X,pConst,scConst);

%% Part 1a/b. Generate Monte Carlo orbit data (nonlinear propagation)
% rng(69420); % Set rng seed for consistency

    % Generate new data by changing the boolean below to true
if false
    fprintf("Generating Monte Carlo Data...\n\n");
    runs = [];
    for k = 1:N
        x0Perturb = chol(P0)*randn(6,1);
    
        data = generateTruthOrbit_MuJ2Drag(pConst, orbital, scConst, x0Perturb, [], tspan);
    
        runs = [runs; data];
    
        fprintf("Generated %.0f orbits so far\n", k)
    end
    
    save(filename, "runs", '-mat');
else
    fprintf("Loading Monte Carlo Data...\n")
    load("..\Data\MonteCarlo3000Runs.mat");
    fprintf("Loaded Data!\n")
end

    % Create nominal orbit and STM
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
fprintf("\nAnalyzing Monte Carlo Data via Nonlinear Propagation...\n")

propTraj_nl = [];
meanTraj_nl = [];
stdTraj_nl = [];
PTraj_nl = {};
for k = 1:length(analysisTimes)
        % Find the index corresponding to the desired analysis time
    idx = find(tspan == analysisTimes(k));

        % Pull out trajectories at this time
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
    titleText = sprintf("Corner Plot at t = %.3f sec, Nonlinear propagation", analysisTimes(k));
    cornerPlot(nom, monteTraj, [], titleText);

        % Save propagated, mean, std, and covariance matrix
    propTraj_nl = [propTraj_nl, nom];
    meanTraj_nl = [meanTraj_nl, mu];
    stdTraj_nl = [stdTraj_nl, sigma];
    PTraj_nl = [PTraj_nl, {P}];

        % Report results
    fprintf("\n\tSummary at t = %.3f, nl prop:\n", analysisTimes(k));
            % Standard Deviations
    fprintf("\t\tState component standard deviations:\n")
    fprintf("\t\t\tX: %.3f,\tY: %.3f,\tZ: %.3f,\tXdot: %.3f,\tYdot: %.3f,\tZdot: %.3f\t\n", stdTraj_nl(:,k))
            % Means
    fprintf("\t\tState component means:\n")
    fprintf("\t\t\tX: %.3f,\tY: %.3f,\tZ: %.3f,\tXdot: %.3f,\tYdot: %.3f,\tZdot: %.3f\t\n", meanTraj_nl(:,k))
            % Nominals
    fprintf("\t\tState component propagated nonlinear values:\n")
    fprintf("\t\t\tX: %.3f,\tY: %.3f,\tZ: %.3f,\tXdot: %.3f,\tYdot: %.3f,\tZdot: %.3f\t\n", propTraj_nl(:,k))
end


%% Part 2. Propagate Uncertainty Using LKF
fprintf("\nAnalyzing Monte Carlo Data via LKF Propagation...\n")

    % Propagate STM
X0 = runs(1).X0 - runs(1).x0Perturb;
XPhi_0 = [X0; reshape(eye(6),36,1)];
[t_LKF, XPhi_LKF] = ode45(@(t,XPhi)STMEOM_MuJ2Drag(t,XPhi,pConst,scConst), analysisTimes, XPhi_0, odeset('AbsTol',1e-12,'RelTol',1e-12));

    % Process STM
Phi_LKF = {};
for k = 1:length(t_LKF)
    Phi = reshape(XPhi_LKF(k,7:end),6,6);
    Phi_LKF = [Phi_LKF; {Phi}];
end

    % Make corner plots and statistics
propTraj_LKF = [];
meanTraj_LKF = [];
stdTraj_LKF = [];
PTraj_LKF = {};
for k = 1:length(analysisTimes)
        % Find the index corresponding to the desired analysis time
    idx = find(tspan == analysisTimes(k));

        % Pull out trajectories at this time
    monteTraj = [];
    for kk = 1:length(runs)
        monteTraj = [monteTraj, runs(kk).X_ref(idx,:)'];
    end
    nom = nominal.X_ref(idx,:)';

        % Calculate mean at this time
    mu = mean(monteTraj,2);

        % Calculate Covariance and standard deviation of trajectories
    P = Phi_LKF{k}*P0*Phi_LKF{k}';
    sigma = [];
    for kk = 1:size(P,1)
        sigma = [sigma; sqrt(P(kk,kk))];
    end

        % Make corner plot
    titleText = sprintf("Corner Plot at t = %.3f sec, LKF propagation", analysisTimes(k));
    cornerPlot(nom, monteTraj, P, titleText);

        % Save propagated, mean, std, and covariance matrix
    propTraj_LKF = [propTraj_LKF, nom];
    meanTraj_LKF = [meanTraj_LKF, mu];
    stdTraj_LKF = [stdTraj_LKF, sigma];
    PTraj_LKF = [PTraj_LKF, {P}];

        % Report results
    fprintf("\n\tSummary at t = %.3f, LKF prop:\n", analysisTimes(k));
            % Standard Deviations
    fprintf("\t\tState component standard deviations:\n")
    fprintf("\t\t\tX: %.3f,\tY: %.3f,\tZ: %.3f,\tXdot: %.3f,\tYdot: %.3f,\tZdot: %.3f\t\n", stdTraj_LKF(:,k))
            % Means
    fprintf("\t\tState component means:\n")
    fprintf("\t\t\tX: %.3f,\tY: %.3f,\tZ: %.3f,\tXdot: %.3f,\tYdot: %.3f,\tZdot: %.3f\t\n", meanTraj_LKF(:,k))
            % Nominals
    fprintf("\t\tState component propagated LKF values:\n")
    fprintf("\t\t\tX: %.3f,\tY: %.3f,\tZ: %.3f,\tXdot: %.3f,\tYdot: %.3f,\tZdot: %.3f\t\n", propTraj_LKF(:,k))
end

%% Part 3. Propagate Uncertainty using UKF
fprintf("\nAnalyzing Monte Carlo Data via UKF Propagation...\n")

    % Compute UKF weights
alpha = 0.5; beta = 2; 
L = length(X0);

kappa = 3 - L;
lambda = (alpha^2)*(L + kappa) - L;
gamma = sqrt(L + lambda);

W_0m = lambda/(L + lambda);
W_0c = lambda/(L + lambda) + (1 - alpha^2 + beta);
W_im = [W_0m, (1/(2*(L + lambda)))*ones(1,2*L)];
W_ic = [W_0c, (1/(2*(L + lambda)))*ones(1,2*L)];
    
    % Propagate sigma points and covariance
X_im1 = X0;
P_im1 = P0;
propTraj_UKF = X0;
PTraj_UKF = P0;
for k = 2:length(analysisTimes)
        % Calculate previous sigma points
    sqrtP_im1 = sqrtm(P_im1); % Used to be chol()
    Chi_im1 = [X_im1, X_im1 + gamma*sqrtP_im1, X_im1 - gamma*sqrtP_im1]; % L x (2L + 1) matrix
    
        % Propagate previous sigma points through dynamics
    ChiVec_im1 = reshape(Chi_im1, L*(2*L+1), 1);
    tspan_UKF = [analysisTimes(k-1), analysisTimes(k)];
    [~,ChiVec] = ode45(@(t,ChiVec)sigPointEOM(t,ChiVec,DynFunc), tspan_UKF, ChiVec_im1, odeset('AbsTol',1e-12,'RelTol',1e-12));
    Chi_i = reshape(ChiVec(end,:), L, 2*L + 1);
    
        % Time update
    X_i = 0;
    for kk = 1:2*L+1
        X_i = X_i + W_im(kk)*Chi_i(:,kk);
    end
    
    P_i = zeros(size(P0));
    for kk = 1:2*L+1
        P_i = P_i + W_ic(kk)*(Chi_i(:,kk) - X_i)*(Chi_i(:,kk) - X_i)';
    end

        % Save values
    propTraj_UKF = [propTraj_UKF, X_i];
    PTraj_UKF = [PTraj_UKF, {P_i}];

        % Update for next time
    X_im1 = X_i;
    P_im1 = P_i;
end

    % Make corner plots and statistics
meanTraj_UKF = [];
stdTraj_UKF = [];
for k = 1:length(analysisTimes)
        % Find the index corresponding to the desired analysis time
    idx = find(tspan == analysisTimes(k));

        % Pull out trajectories at this time
    monteTraj = [];
    for kk = 1:length(runs)
        monteTraj = [monteTraj, runs(kk).X_ref(idx,:)'];
    end

        % Calculate mean at this time
    mu = mean(monteTraj,2);

        % Pull out covariance and propagated trajectory at this time
    nom = propTraj_UKF(:,k);
    P = PTraj_UKF{k};

        % Calculate standard deviations
    sigma = [];
    for kk = 1:size(P,1)
        sigma = [sigma; sqrt(P(kk,kk))];
    end

        % Make corner plot
    titleText = sprintf("Corner Plot at t = %.3f sec, UKF propagation", analysisTimes(k));
    cornerPlot(nom, monteTraj, P, titleText);

        % Save nominal, mean, std, and covariance matrix
    meanTraj_UKF = [meanTraj_UKF, mu];
    stdTraj_UKF = [stdTraj_UKF, sigma];

        % Report results
    fprintf("\n\tSummary at t = %.3f, UKF prop:\n", analysisTimes(k));
            % Standard Deviations
    fprintf("\t\tState component standard deviations:\n")
    fprintf("\t\t\tX: %.3f,\tY: %.3f,\tZ: %.3f,\tXdot: %.3f,\tYdot: %.3f,\tZdot: %.3f\t\n", stdTraj_UKF(:,k))
            % Means
    fprintf("\t\tState component means:\n")
    fprintf("\t\t\tX: %.3f,\tY: %.3f,\tZ: %.3f,\tXdot: %.3f,\tYdot: %.3f,\tZdot: %.3f\t\n", meanTraj_UKF(:,k))
            % Nominals
    fprintf("\t\tState component propagated UKF values:\n")
    fprintf("\t\t\tX: %.3f,\tY: %.3f,\tZ: %.3f,\tXdot: %.3f,\tYdot: %.3f,\tZdot: %.3f\t\n", propTraj_UKF(:,k))
end

%% Part 4. Propagate Uncertainty using Gaussian Sums
fprintf("\nAnalyzing Monte Carlo Data via Gaussian Sums...\n")







