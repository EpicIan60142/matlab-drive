%% ECEN 4138 HW 1 Problem 4 Script
%   By: Ian Faber, 09/06/2023

%% Housekeeping
clc; clear; close all;

%% Setup
m = 1500; % kg
b = 70; % N-sec/m
K = 20000; % N-sec/m, vary for this problem

vr = 1; % m/s

tspan = 0:0.001:10;

titleText = sprintf("Car velocity vs. time, K = %.0f N/m/s", K);

disturbTime = tspan(end)/2;

const = [m; b; K; vr; disturbTime];

v0 = 0; % IC vector

%% Simulate

[time, state] = ode45(@(t,state)cruiseEOM(t,state,const), tspan, v0);

%% Analyze

v = state(:,1);

idx = find(time == disturbTime);
vrPlot = zeros(size(time));
vrPlot(idx:end) = vr;

figure
hold on; grid on;
title(titleText)
plot(time, v, 'b-')
plot(time, vrPlot, 'k--')
xlabel("Time [sec]")
ylabel("v [m/s]")

legend("Car velocity", "Velocity command", 'Location', 'best')

%% EOM function
function dX = cruiseEOM(t,X,const)
% EOM function for simulating a cruise control system with ode45
%   Inputs:
%       t: time [sec]
%       X: state vector
%           v
%       const: vector of constants for simulation
%           [m; b; K; vr; disturbTime]
%
%   Outputs:
%       dX: rate of change vector
%           [ vx; vy; ax; ay ]
%
%   By: Ian Faber, 09/05/2023
%

    v = X;

    m = const(1);
    b = const(2);
    K = const(3);
    vr = const(4);
    disturbTime = const(5);

    if t < disturbTime
        vr = 0;
    end

    a = (1/m)*(K*vr - (b+K)*v);

    dX = a;

end
