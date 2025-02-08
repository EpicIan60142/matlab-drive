%% Housekeeping
clc; clear; close all;

%% Constants
ratio = 2.5/(3.959e-5);

lb2n = 4.44882; % lbf to n
in2m = 0.0254; % in to m

%% Data Extraction
trussData = readmatrix("Data/Test 02");
loading = lb2n*trussData(:,1);
F0 = lb2n*trussData(:,2);
F1 = lb2n*trussData(:,3);
F2 = lb2n*trussData(:,4);
F3D = lb2n*trussData(:,5);
LVDT = in2m*trussData(:,6);

ReactionForces = load("Data\ReactionForces1.txt");
BarForces = load("Data\BarForces1.txt");

%% Regression & Uncertainty
% F0 vs. loading
[coef0, approxCurve0] = leastSquares(loading, F0, 1);

% F1 vs. loading
[coef1, approxCurve1] = leastSquares(loading, F1, 1);

% F2 vs. loading
[coef2, approxCurve2] = leastSquares(loading, F2, 1);

% F3D vs. loading
[coef3D, approxCurve3D] = leastSquares(loading, F3D, 1);

% LVDT vs. loading
[coefLV, approxCurveLV] = leastSquares(loading, LVDT, 1);

% F0 vs. curve
diffF0 = F0 - approxCurve0(loading);
uncertaintyF0 = [mean(diffF0), std(diffF0)];

% F1 vs. curve
diffF1 = F1 - approxCurve1(loading);
uncertaintyF1 = [mean(diffF1), std(diffF1)];

% F2 vs. curve
diffF2 = F2 - approxCurve2(loading);
uncertaintyF2 = [mean(diffF2), std(diffF2)];

% F3D vs. curve
diffF3D = F3D - approxCurve3D(loading);
uncertaintyF3D = [mean(diffF3D), std(diffF3D)];

% LVDT vs. curve
diffLVDT = LVDT - approxCurveLV(loading);
uncertaintyLVDT = [mean(diffLVDT), std(diffLVDT)];

%% FEM Setup
% The order for these are "node","force value"
ReactF = [1, 55.562; 52, 55.638; 17, 55.638; 68, 55.562];
BarF = [128, -440.02; 131, -440.02; 178, -388.15; 180, -388.22];

% The order for this is "node","displacement"
NodeDisp = [26,ratio*(-0.29207e-007);43,ratio*(-0.29263e-007)];

% Info for Figures
BarFfront = linspace(0,1,10)*351.233 - 400.34; % force in bar on bottom front of truss
dispMiddle = linspace(0,1,10)*(ratio*(-0.29207e-007)); % Displacement of nodes in the middle of truss  
extLoad = linspace(0,1,10)*(222.4); % external load value (total)
ReactF1 = linspace(0,1,10)*55.562; % reaction force at node 1 (pin)
ReactF52 = linspace(0,1,10)*55.638; % reaction force at node 52 (pin)
ReactFRoller = linspace(0,1,10)*(55.638+55.562);

%% Comparison
figure(1)
hold on
grid on
title("Reaction Force F0, Comparison Between FEM and Experiment")
plot(loading, approxCurve0(loading), 'b-')
plot(extLoad, ReactF1, 'r-')
xlabel("External Loading (N)")
ylabel("Reaction Force F0 (N)")
legend("Experiment", "FEM", 'Location', 'best')

figure(2)
hold on
grid on
title("Reaction Force F1, Comparison Between FEM and Experiment")
plot(loading, approxCurve1(loading), 'b-')
plot(extLoad, ReactF52, 'r-')
xlabel("External Loading (N)")
ylabel("Reaction Force F1 (N)")
legend("Experiment", "FEM", 'Location', 'best')

figure(3)
hold on
grid on
title("Reaction Force F2, Comparison Between FEM and Experiment")
plot(loading, approxCurve2(loading), 'b-')
plot(extLoad, ReactFRoller, 'r-')
xlabel("External Loading (N)")
ylabel("Reaction Force F2 (N)")
legend("Experiment", "FEM", 'Location', 'best')

figure(4)
hold on
grid on
title("Internal Force F3D, Comparison Between FEM and Experiment")
plot(loading, approxCurve3D(loading), 'b-')
plot(extLoad, BarFfront, 'r-')
xlabel("External Loading (N)")
ylabel("Internal Force F3D (N)")
legend("Experiment", "FEM", 'Location', 'best')

figure(5)
hold on
grid on
title("Vertical Displacement, Comparison Between FEM and Experiment")
plot(loading, approxCurveLV(loading), 'b-')
plot(extLoad, dispMiddle, 'r-')
xlabel("External Loading (N)")
ylabel("Vertical Displacement (m)")
legend("Experiment", "FEM", 'Location', 'best')





