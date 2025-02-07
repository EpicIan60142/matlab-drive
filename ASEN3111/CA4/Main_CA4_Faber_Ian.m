%% ASEN 3111 - Computational Assignment 04 - Main
%   Driver program for CA 04 - Compressible Aerodynamics
%       Program that models compressible aerodynamics for thin diamond
%       airfoils and compares Shock Expansion Theory with Linear Supersonic
%       Flow. Also verifies the Theta-Beta-Mach relation by recreating
%       Anderson Figure 9.9.
%
%   Author: Ian Faber
%   Collaborators: Maggie Wussow
%   Date started: 04/13/2023
%   Date finished: 04/26/2023

%% Housekeeping
clc; clear; close all;

tic; % Start timing

%% Starter Command Window Management
fprintf("\t+--------------------------------+\n")
fprintf("\t\t CA 4 Results - Ian Faber\n\n")

%% Problem 1
% fprintf("\n\t Problem 1:\n\n")

% Constants
gamma = 1.4; % air

small = 50; % Strong shock beta cutoff for plotting

startCol = [1,0,0]; % Color for small mach numbers
stopCol = [0,0,1]; % Color for large mach numbers

% Setup
maxTheta = 45.5; % degrees

nTheta = 5000; % Number of theta points to solve at

theta = linspace(0, maxTheta, nTheta); % Define theta vector

M = [(2.2:0.2:4.0)'; 4.5; 5; (6:2:10)'; 20]; % Define mach vector

% Preallocate beta outputs for speed
%   - Organized in a 3-d matrix as follows: ([theta; beta] x [num thetas] x
%     [num machs]
%   - For every mach, calculate 2 vectors (theta and beta) that are as long
%     as the specified theta test vector
betaWeak = zeros(2,nTheta,length(M));
betaStrong = zeros(2,nTheta,length(M));

% Calculate betas and tabulate them if they satisfy physical criteria
for k = 1:length(M)
    for kk = 1:length(theta)
        bWeak = ObliqueShockBeta(M(k), theta(kk), gamma, 'Weak');
        bStrong = ObliqueShockBeta(M(k), theta(kk), gamma, 'Strong');

        if isreal(bWeak) && ~isnan(bWeak) && bWeak > 0 
            betaWeak(1,kk,k) = theta(kk);
            betaWeak(2,kk,k) = bWeak;
        end

        if isreal(bStrong) && ~isnan(bStrong) && bStrong > 0
            betaStrong(1,kk,k) = theta(kk);
            betaStrong(2,kk,k) = bStrong;
        end

    end
end

% Recreate Fig. 9.9 from Anderson (AKA NACA-TR-1135)
figure
hold on
grid on
title("$\theta$ - $\beta$ - M Relation",'Interpreter','latex')

for k = 1:length(M)
    color = map(k, 1, length(M), startCol, stopCol);

    stop = find(betaStrong(2,:,k) > small, 1, 'last');

    curve(k) = plot(betaWeak(1,1:stop,k), betaWeak(2,1:stop,k), 'Color', color); % Plot strong beta (y) vs. theta (x) for a given mach number
    plot(betaStrong(1,1:stop,k),betaStrong(2,1:stop,k), 'Color', color); % Plot weak beta (y) vs. theta (x) for a given mach number
    
    label(k) = sprintf("M = %.1f", M(k));
end

xlabel("Wedge Angle \theta [deg]")
ylabel("Shock Angle \beta [deg]")
xlim([26,58])
ylim([0,90])

legend(curve, label)

%% Problem 2
fprintf("\n\tProblem 2:\n\n")

% Geometry to investigate
M = 3;
alpha = 10; % deg
epsilon1 = 7.5; % deg
epsilon2 = 5; % deg

[c_l, c_dw] = DiamondAirfoil(M, alpha, epsilon1, epsilon2);

fprintf("Sectional lift coefficient: c_l = %.4f\n", c_l)
fprintf("Sectional wave drag coefficient: c_dw = %.4f\n", c_dw)


%% Problem 3

fprintf("\n\tProblem 3:\n\n")

% Define study space
alpha = -20:0.1:20;
M = [2, 3, 4, 5];
color = ['b','r','g','m'];

% Define gl^2 and gu^2
l1c = 1/(cosd(epsilon1) + (sind(epsilon1)/tand(epsilon2))); % L1/c
x1c = l1c*cosd(epsilon1); % x1/c, location where l1 and l2 intersect

gl2 = x1c*(tand(epsilon1)^2 - tand(epsilon2)^2) + tand(epsilon2)^2;

gu2 = x1c*(tand(epsilon1)^2 - tand(epsilon2)^2) + tand(epsilon2)^2; % Same as gl2 because of vertical symmetry

% Allocate matrices for storing results
c_lStudy = zeros(length(alpha), length(M));
c_dwStudy = zeros(length(alpha), length(M));
c_lLinear = zeros(length(alpha), length(M));
c_dwLinear = zeros(length(alpha), length(M));

% Calculate c_l and c_dw for each test case and their linear counterparts
for k = 1:length(M)
    for kk = 1:length(alpha)
        [c_lStudy(kk,k), c_dwStudy(kk,k)] = DiamondAirfoil(M(k), alpha(kk), epsilon1, epsilon2);
        c_lLinear(kk,k) = (4*deg2rad(alpha(kk)))/sqrt(M(k)^2 - 1);
        c_dwLinear(kk,k) = (2/sqrt(M(k)^2 - 1))*(2*deg2rad(alpha(kk))^2 + gl2 + gu2);
    end
end

% Plot c_l
figure
hold on
grid on
title("c_l vs. \alpha at Various Mach Numbers")
for k = 1:length(M)
    machs(k) = plot(alpha, c_lStudy(:,k), join([color(k),'-'],''));
    linears(k) = plot(alpha, c_lLinear(:,k), join([color(k),'--'],''));
    machLabel(k) = sprintf("S.E.T., M = %.0f", M(k));
    linearLabel(k) = sprintf("L.S.F., M = %.0f", M(k));
end
xlabel("\alpha [deg]")
ylabel("c_l")
subset = [machs, linears];
label = [machLabel, linearLabel];
legend(subset, label, 'Location','best', 'NumColumns', 2)

% Plot c_dw
figure
hold on
grid on
title("c_{d,w} vs. \alpha at Various Mach Numbers")
for k = 1:length(M)
    machs(k) = plot(alpha, c_dwStudy(:,k), join([color(k),'-'],''));
    linears(k) = plot(alpha, c_dwLinear(:,k), join([color(k),'--'],''));
    machLabel(k) = sprintf("S.E.T., M = %.0f", M(k));
    linearLabel(k) = sprintf("L.S.F., M = %.0f", M(k));
end
xlabel("\alpha [deg]")
ylabel("c_{d,w}")
subset = [machs, linears];
label = [machLabel, linearLabel];
legend(subset, label, 'Location','best', 'NumColumns', 2)

%% End Code
fprintf("\n")
toc; % Stop timing

%% Final Command Window Management
fprintf("\n\n\t+--------------------------------+\n\n")
disp(getRalphie)


