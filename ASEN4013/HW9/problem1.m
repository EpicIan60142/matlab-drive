%% ASEN 4013 HW 9 Problem 1
% - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
P = (1:1:10)*10^6; % MPa
Pc = P.^0.55; % MPa

g0 = 9.807; % m/s^2

%% Data
cStar = [1612.4, 1622.2, 1627.5, 1631.1, 1633.7, 1635.8, 1637.5, 1638.9, 1640.1, 1641.2]; % m/s

Isp = [208.65, 233.68, 245.93, 253.72, 259.31, 263.61, 267.06, 269.94, 272.38, 274.50]; % sec

cStarAvg = mean(cStar)

IspAvg = mean(Isp)

%% Plotting
figure
hold on
plot(P, cStar, '.', 'MarkerSize', 20);
% plot(P, mean(coefCStar).*Pc, 'k--');

figure
hold on
plot(P, Isp, '.', 'MarkerSize', 20);
% plot(P, coefIsp*Isp, 'k--');




