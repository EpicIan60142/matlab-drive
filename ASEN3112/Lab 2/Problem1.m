%% ASEN 3112 Lab 2 Problem 1 Script
% Section 013 Group 2 - Ian Faber, Ada Forsner, Ethan Breaux, Matt Bechtel, Rithik Gangopadhyay

%% Housekeeping
clc; clear; close all;

%% Constants
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

%% Regression
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

%% Uncertainty
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

%% Plotting
figure
hold on
grid on
title("F0 vs. External Loading")
plot(loading, F0, 'b.')
plot(loading, approxCurve0(loading), 'k--')
xlabel("External Loading (N)")
ylabel("F0 (N)")
legend("Raw Data", "Best Fit")

figure
hold on
grid on
title("F1 vs. External Loading")
plot(loading, F1, 'b.')
plot(loading, approxCurve1(loading), 'k--')
xlabel("External Loading (N)")
ylabel("F1 (N)")
legend("Raw Data", "Best Fit")

figure
hold on
grid on
title("F2 vs. External Loading")
plot(loading, F2, 'b.')
plot(loading, approxCurve2(loading), 'k--')
xlabel("External Loading (N)")
ylabel("F2 (N)")
legend("Raw Data", "Best Fit")

figure
hold on
grid on
title("F3D vs. External Loading")
plot(loading, F3D, 'b.')
plot(loading, approxCurve3D(loading), 'k--')
xlabel("External Loading (N)")
ylabel("F3D (N)")
legend("Raw Data", "Best Fit")

figure
hold on
grid on
title("LVDT vs. External Loading")
plot(loading, LVDT, 'b.')
plot(loading, approxCurveLV(loading), 'k--')
xlabel("External Loading (N)")
ylabel("LVDT (m)")
legend("Raw Data", "Best Fit")



