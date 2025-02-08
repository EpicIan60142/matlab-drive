%% ECEN 3300 HW 3 Bonus Question
%       by: Ian Faber

% Housekeeping
clc; clear; close all;

% Code
n = -20:1:20;
x = @(n) sin((pi/8)*n);

x1 = x(n);
y1 = x(2*n);
y2 = x(n/2);

figure
hold on
title("x[n] manipulations")
plot(n,x1);
plot(n,y1);
plot(n,y2);
ylim([-2,2]);
legend("x[n]", "x[2n]", "x[n/2]", 'location','best')