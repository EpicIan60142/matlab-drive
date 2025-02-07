% Ian Faber, Ashton Miner, Teegan Oatley, Chaney Sullivan
% ASEN 3128-011
% ASEN3128Lab1Problem1.m
% Created: 8/23/22

%% Housekeeping
clc; clear; close all;

%% ODE45 Setup

% Initial conditions
X0 = [1; 1; 0.01];

% Integration time span
tspan = [0 .7];

% Tolerance settings
setup = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

%% ODE45 Call
[t, state] = ode45(@(t, state)EOM(t, state), tspan, X0, setup);

%% Plotting
figure(1)
sgtitle("1.a. ODE45 Integration Output")
subplot(3,1,1)
hold on
title("X vs Time")
xlabel("Time")
ylabel("X")
plot(t, state(:,1))
hold off

subplot(3,1,2)
hold on
title("Y vs Time")
xlabel("Time")
ylabel("Y")
plot(t, state(:,2))
hold off

subplot(3,1,3)
hold on
title("Z vs Time")
xlabel("Time")
ylabel("Z")
plot(t, state(:,3))
hold off
