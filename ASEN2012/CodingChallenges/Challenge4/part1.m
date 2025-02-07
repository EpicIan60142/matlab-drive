%% Ian Faber - SID: 108577813
% ASEN 2012: Coding Challenge 4
% Last modified: 10/14/21


%% Housekeeping

clc;         % Clear command window
close all;   % Close out of any open figures
hold on;

%% Define the Analytical Solution

% The expression of the exact solution
f_definite = (1-exp(-5).*(sin(5)+cos(5)))/2

%% Step 1
% Import 'DataPoints.mat' 
data = load('DataPoints.mat');
x_1 = data.x_1;
x_2 = data.x_2;
x_3 = data.x_3;
y_1 = data.y_1;
y_2 = data.y_2;
y_3 = data.y_3;

%% Step 2
% Generate Riemann sums with the coarsest spacing. Note that it is possible to calculate 
% these sums in a single line, though it is possible with for-loops as well.
dx = x_1(2);
soln_left_1 = sum(y_1(1:(length(y_1)-1))*dx)
soln_right_1 = sum(y_1(1:length(y_1))*dx)

% Generate other solutions the same way with the other two resolutions. Ideally, these should be 
% identical as the above, except for the x and y data used.
dx = x_2(2);
soln_left_2 = sum(y_2(1:(length(y_2)-1))*dx)
soln_right_2 = sum(y_2(1:length(y_2))*dx)

dx = x_3(2);
soln_left_3 = sum(y_3(1:(length(y_3)-1))*dx)
soln_right_3 = sum(y_3(1:length(y_3))*dx)

plot(x_3, y_3, 'k');


%% Step 3
% Generate your own data points from a new function to estimate the integral.
n = 100000
x_4 = linspace(0,1,n);
y_4 = sqrt(1-(x_4.^4));

dx = 1/n;
soln_4 = sum(y_4(1:length(y_4))*dx)

% Leave a comment if you can identify if this is an over or under estimate
%plot(x_4, y_4, 'k')
% This is an underestimate, as the function is concave down throughout the interval of interest
