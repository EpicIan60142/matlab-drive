%% ASEN 5010 HW 3 Concept Check 8 Script
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Problem
beta0 = [0.408248; 0; 0.408248; 0.816497];

t = (0:0.001:60)';

[time, beta] = ode45(@(t, state)betaEOM(t, state), t, beta0);

timeInstant = 42;

answer = norm(beta(time == timeInstant,2:4))

figure
hold on
sgtitle("Euler Parameter Evolution")

subplot(4,1,1)
plot(time, beta(:,1));
xlabel("time [sec]")
ylabel("\beta_0")

subplot(4,1,2)
plot(time, beta(:,2));
xlabel("time [sec]")
ylabel("\beta_1")

subplot(4,1,3)
plot(time, beta(:,3));
xlabel("time [sec]")
ylabel("\beta_2")

subplot(4,1,4)
plot(time, beta(:,4));
xlabel("time [sec]")
ylabel("\beta_3")


