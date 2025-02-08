%% ASEN 5010 HW 1 CC 9 Q2 Script
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
theta0 = [40; 30; 80]; % deg

t = (0:0.001:60)'; 

[time, angle] = ode45(@(t,state)angleEOM(t,state), t, theta0);

answer = norm(deg2rad(angle((time==42),:)))
 


