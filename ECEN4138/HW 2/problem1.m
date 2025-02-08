%% ECEN 4138 HW 2 Problem 1 Script
%   By: Ian Faber, 09/11/2023

%% Housekeeping
clc; clear; close all;

%% Setup
M = 1; % kg
m = 1; % kg
b1 = 0.1; % N/m/s
b2 = [0, 0.1, 1]; % N/m/s, to turn off air drag b2 = 0
k = 1; % N/m

tspan = 0:0.001:50; % sec

u = 1; % N, input force

x0 = [0; 0; 0; 0]; % x; y; vx; vy

%% Simulate

for ii = 1:length(b2)
    

%     titleTextX = sprintf("Flexible system x response with b_2 = %.3f", b2(ii));
%     titleTextY = sprintf("Flexible system y response with b_2 = %.3f", b2(ii));
%     titleTextVx = sprintf("Flexible system v_x response with b_2 = %.3f", b2(ii));
    titleTextVy = sprintf("Flexible system v_y response with b_2 = %.3f", b2(ii));

    const = [M; m; b1; b2(ii); k; u];

    [time, state] = ode45(@(time, state)flexEOM(time, state, const), tspan, x0);

    %% Analyze
    x = state(:,1);
    y = state(:,2);
    vx = state(:,3);
    vy = state(:,4);
    
%     figure
%     hold on; grid on;
%     title(titleTextX)
%     plot(time, x, 'b-')
%     xlabel("Time [sec]")
%     ylabel("x position [m]")
% 
%     figure
%     hold on; grid on;
%     title(titleTextY)
%     plot(time, y, 'b-')
%     xlabel("Time [sec]")
%     ylabel("y position [m]")
% 
%     figure
%     hold on; grid on;
%     title(titleTextVx)
%     plot(time, vx, 'b-')
%     xlabel("Time [sec]")
%     ylabel("x velocity [m/s]")

    figure
    hold on; grid on;
    title(titleTextVy)
    plot(time, vy, 'b-')
    xlabel("Time [sec]")
    ylabel("y velocity [m/s]")

end

%% EOM function

function dX = flexEOM(t, X, const)
% EOM function for simulating a simplified 1-D flexible system with ode45
%   Inputs:
%       t: time [sec]
%       X: state vector
%           [ x; y; vx; vy ]
%       const: vector of constants for simulation
%           [M; m; b1; b2; k; u] -> If not simulating air drag, b2 = 0
%
%   Outputs:
%       dX: rate of change vector
%           [ vx; vy; ax; ay ]
%
%   By: Ian Faber, 09/11/2023
%

    M = const(1);
    m = const(2);
    b1 = const(3);
    b2 = const(4);
    k = const(5);
    u = const(6);
    
    x = X(1);
    y = X(2);
    vx = X(3);
    vy = X(4);
    
%     if t > 30
%         u = 0;
%     end
        
    ax = (1/m)*(b1*vy + k*y - b1*vx - k*x);
    ay = (1/M)*(u + b1*vx + k*x - (b1+b2)*vy - k*y );
    
    dX = [vx; vy; ax; ay];

end
