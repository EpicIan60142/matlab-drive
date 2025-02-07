clc; clear; close all;

t = 0:0.0001:2;

u1 = @(t) (12/5)*cos(((10*sqrt(15))/3)*t) - (2/5)*cos(10*sqrt(10)*t);
u2 = @(t) (8/5)*cos(((10*sqrt(15))/3)*t) + (2/5)*cos(10*sqrt(10)*t);
% uMean = @(t) (u1(t) + u2(t))/2;

figure
hold on
title("Problem 9.3. u_1(t), u_2(t) vs. time")
plot(t, u1(t), 'b-')
plot(t, u2(t), 'r-')
% plot(t, uMean(t), 'k-')
xlabel("Time (sec)")
ylabel("Response (in)")

legend("u_1(t)", "u_2(t)", 'Location', 'best')





