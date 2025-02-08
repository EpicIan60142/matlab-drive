%% ASEN 3113 Prelab 2 Question 3: Experimental Steady State Solution
%   - Ian Faber

%% Housekeeping
clc; clear; close all

%% Code
Th = [17.6, 21.61, 25.13, 29.22, 34.92, 38.10, 45.21, 47.01];
x = 1.375:0.5:1.375+7*0.5;

xLine = 0:0.5:x(end)+0.5;

coef = polyfit(x,Th,1);

T = polyval(coef, xLine); 

T0 = polyval(coef,0)
H = coef(1,1)

TLine = polyval(coef,x);
error = 100*(Th - TLine)./Th;

lineLabel = sprintf("Line of best fit: T = %.3f*x + %.3f", coef(1,1), coef(1,2))

figure
hold on
title("Temperature vs. Linear Distance")
plot(x,Th, '.', 'MarkerSize',15)
plot(xLine,T,'k--')
xlabel("Distance [in]")
ylabel("Temperature [^oC]")
legend("Thermocouple data",lineLabel,'Location','best')

figure
hold on
title("Temperature Best Fit Error")
stem(x,error)
xlabel("Distance [in]")
ylabel("Percent Error [%]")
