%% ECEN 3300 HW 4 Problem 6 Plotting
% Ian Faber

%% Housekeeping
clc; clear; close all

%% Code
n = 0:1:100;

func = @(n) (n>=0).*0.5.^n + ((n-1)>=0).*0.5.^(n-1) + ((n-2)>=0).*0.5.^(n-2) - 0.5.^(n-5) + 4;

figure
hold on
title("ECEN 3300 HW 4 Problem 6 y[n] Plot")
plot(n, func(n));
xlabel("n")
ylabel("y[n]")