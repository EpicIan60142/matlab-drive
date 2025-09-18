%% ASEN 5254 HW 3 Main Script
% By: Ian Faber, 09/17/2025

%% Housekeeping
clc; clear; close all;

%% Setup
testC3 = -1:0.01:1; % List of cos(theta3) to test
point = [0; 4]; % 

%% Run procedure to solve for thetas and check answer
thetas = solveThetas(point, testC3);

for k = 1:length(thetas.theta1)
    theta1 = thetas.theta1(k);
    theta2 = thetas.theta2(k);
    theta3 = thetas.theta3(k);

    x = cos(theta1)*(cos(theta2)*(9*cos(theta3) + 8) - 9*sin(theta2)*sin(theta3) + 8) ...
      + sin(theta1)*(sin(theta2)*(-9*cos(theta3) - 8) - 9*cos(theta2)*sin(theta3));

    y = cos(theta1)*(sin(theta2)*(9*cos(theta3) + 8) + 9*cos(theta2)*sin(theta3)) ...
      + sin(theta1)*(cos(theta2)*(9*cos(theta3) + 8) - 9*sin(theta2)*sin(theta3) + 8);

    fprintf("End effector with theta = [%.3f, %.3f, %.3f]: (x,y) = [%.3f, %.3f]\n", theta1, theta2, theta3, x, y)
end

%% Plot results
fig = figure;
hold on; grid on;
title("Valid solutions for \theta_1, \theta_2, \theta_3")
plot3(rad2deg(thetas.theta1), rad2deg(thetas.theta2), rad2deg(thetas.theta3), '.');
xlabel("\theta_1 [deg]"); ylabel("\theta_2 [deg]"); zlabel("\theta_3 [deg]");
view([30 35]);
