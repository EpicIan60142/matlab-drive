%% ASEN 5254 HW 8 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
nAgents = [2, 3, 4, 5, 6];

centralTimes = [996.117, 3793.81, 14954.6, 16441.7, 15623.6]; % ms
centralSizes = [569.13, 1184.87, 2968.75, 2305.71, 1537.78];

decoupTimes = [611.104, 915.535, 1514.47, 1826.3, 2247.91]; % ms

markerSize = 20;

%% Plot results

figure;
hold on; grid on;
title("Centralized Average Compute Time vs. Number of Agents")
plot(nAgents, centralTimes, 'b.', 'MarkerSize', markerSize);
xlabel("Number of Agents"); ylabel("Average compute time [ms]")

figure;
hold on; grid on;
title("Centralized Average Tree Size vs. Number of Agents")
plot(nAgents, centralSizes, 'r.', 'MarkerSize', markerSize);
xlabel("Number of Agents"); ylabel("Average tree size")

figure;
hold on; grid on;
title("Decoupled Average Compute Time vs. Number of Agents")
plot(nAgents, decoupTimes, 'b.', 'MarkerSize', markerSize);
xlabel("Number of Agents"); ylabel("Average compute time [ms]")

figure;
hold on; grid on;
title("Average Compute Time Comparison vs. Number of Agents")
plot(nAgents, centralTimes, 'b.', 'MarkerSize', markerSize);
plot(nAgents, decoupTimes, 'r.', 'MarkerSize', markerSize);
xlabel("Number of Agents"); ylabel("Average compute time [ms]")
legend("Centralized", "Decoupled", 'location', 'best');
