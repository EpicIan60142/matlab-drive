%% ECEN 4138 HW 1 Problem 1 Script
%   By: Ian Faber, 09/05/2023

%% Housekeeping
clc; clear; close all;

%% Setup
m = 1; % kg
L = 1; % m
g = 9.8; % m/s^2

tspan = [0, 10]; % Simulate from 0 to 10 seconds

theta0 = deg2rad([1,5,25]); % deg -> rad
omega0 = deg2rad([0,0,0]); % deg/s -> rad/s

timeL = cell(1,length(theta0));
stateL = cell(1,length(theta0));
timeNL = cell(1,length(theta0));
stateNL = cell(1,length(theta0));

chartText = strings(length(theta0),1);
thetaText = strings(length(theta0),1);
omegaText = strings(length(theta0),1);

x0 = [theta0; omega0];

%% Simulate

for k = 1:length(theta0)

    x0 = [theta0(k); omega0(k)];

    chartText(k) = sprintf("Linear vs. Nonlinear Pendulum Simulation \n \\theta_0 = %.3f^o and \\omega_0 = %.3f deg/s", rad2deg(theta0(k)), rad2deg(omega0(k)));
    thetaText(k) = sprintf("\\theta vs. time");
    omegaText(k) = sprintf("\\omega vs. time");

    [tL, sL] = ode45(@(t,state)pendulumEOM(t,state,g,0), tspan, x0); % Linear EOM

    [tNL, sNL] = ode45(@(t,state)pendulumEOM(t,state,g,1), tspan, x0); % Nonlinear EOM

    timeL{k} = tL;
    stateL{k} = sL;

    timeNL{k} = tNL;
    stateNL{k} = sNL;

end

clear tL sL tNL sNL

%% Analyze

for k = 1:length(theta0)

    tL = timeL{:,k};
    sL = stateL{:,k};

    tNL = timeNL{:,k};
    sNL = stateNL{:,k};

    thetaL = sL(:,1);
    omegaL = sL(:,2);

    thetaNL = sNL(:,1);
    omegaNL = sNL(:,2);
    

    figure
    sgtitle(chartText(k))

    subplot(2,1,1)
    hold on; grid on;
    title(thetaText(k))
    plot(tL, rad2deg(thetaL), 'b-')
    plot(tNL, rad2deg(thetaNL), 'r-')
    xlabel("Time [sec]")
    ylabel("\theta [deg]")
    
    legend("Linear", "Nonlinear")

    subplot(2,1,2)
    hold on; grid on;
    title(omegaText(k))
    plot(tL, rad2deg(omegaL), 'b-')
    plot(tNL, rad2deg(omegaNL), 'r-')
    xlabel("Time [sec]")
    ylabel("\omega [deg/s]")
    
    legend("Linear", "Nonlinear")

end

%% EOM function
function dX = pendulumEOM(t,X,g,config) 
% EOM function for simulating a simple pendulum with ode45
%   Inputs:
%       t: time [sec]
%       X: state vector
%           [ theta; omega ]
%       config: Type of EOM to simulate
%           0 = Linear, 1 = Nonlinear, defaults to linear if not 0 or 1
%
%   Outputs:
%       dX: rate of change vector
%           [ omega; alpha ]
%
%   By: Ian Faber, 09/05/2023
%

% Extract state variables
theta = X(1);
omega = X(2);

% Choose equation set to simulate
switch config
    case 0 % Linear
        alpha = -g*theta;
    case 1 % Nonlinear
        alpha = -g*sin(theta);
    otherwise % Dumb user
        alpha = -g*theta;
end

dX = [omega; alpha];

end
