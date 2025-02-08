%% ECEN 4138 HW 1 Problem 2 Script
%   By: Ian Faber, 09/05/2023

%% Housekeeping
clc; clear; close all;

%% Setup
m1 = 10; % kg
m2 = 350; % kg
Kw = 500000; % N/m
Ks = 10000; % N/m
b = 7500; % N/m/s, change for problem, start at 0

r = 1; % m, unit road disturbance applied at some time during simulation

titleText = sprintf("Car motion vs. time, b = %.0f N/m/s", b);

tspan = 0:0.001:10;

disturbTime = tspan(end)/2; % Time the disturbance is applied in the simulation

const = [m1; m2; Kw; Ks; b; r; disturbTime];

x0 = [0; 0; 0; 0]; % IC vector, formatted [x, y, vx, vy]

%% Simulate

[time, state] = ode45(@(t,state)suspensionEOM(t,state,const), tspan, x0);

%% Analyze

x = state(:,1);
y = state(:,2);
vx = state(:,3);
vy = state(:,4);

idx = find(time == disturbTime);
rPlot = zeros(size(time));
rPlot(idx:end) = r;

figure
hold on; grid on;
title(titleText)
plot(time, y, 'b-')
plot(time, rPlot, 'k--')
xlabel("Time [sec]")
ylabel("y [m]")

legend("Car motion", "Road disturbance", 'Location', 'best')

%% EOM function
function dX = suspensionEOM(t,X,const)
% EOM function for simulating a simplified car suspension with ode45
%   Inputs:
%       t: time [sec]
%       X: state vector
%           [ x; y; vx; vy ]
%       const: vector of constants for simulation
%           [m1; m2; Kw; Ks; b; r; disturbTime]
%
%   Outputs:
%       dX: rate of change vector
%           [ vx; vy; ax; ay ]
%
%   By: Ian Faber, 09/05/2023
%

    x = X(1);
    y = X(2);
    vx = X(3);
    vy = X(4);
    
    m1 = const(1);
    m2 = const(2);
    Kw = const(3);
    Ks = const(4);
    b = const(5);
    r = const(6);
    disturbTime = const(7);
    
    if t < disturbTime % Don't apply the disturbance r until specified in the main script
        r = 0;
    end
    
    ax = (Kw/m1)*(r-x) + (Ks/m1)*(y-x) + (b/m1)*(vy-vx);
    ay = (Ks/m2)*(x-y) + (b/m2)*(vx-vy);
    
    dX = [vx; vy; ax; ay];

end

