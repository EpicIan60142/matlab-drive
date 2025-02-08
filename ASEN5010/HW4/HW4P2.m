%% ASEN 5010 HW 4 Problem 2 Script
%   - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath("..\Utilities\")

sig0 = [0; 0; 0];
w0 = {[1; 0; 0], [1; 0.1; 0], [0; 1; 0], [0; 1; 0.1], [0; 0; 1], [0; 0.1; 1]}; % Simulation cases

I = [   125,    0,      0;
        0,      100,    0;
        0,      0,      75
    ];
    
dt = 0.01; % time step
t = (0:dt:60)'; % Simulate for 1 minute

% Starter variables for plotting
wMax = 0;
wMin = 9999999;
buffer = 1.25; % Plotting buffer

%% Run simulations

for k = 1:length(w0)
    x0 = [sig0; w0{k}]; % Update simulation case
    u0 = zeros(3,1);
    
    % Run RK4 algorithm
    output{k} = RK4_RigidBody_MRP(x0, u0, I, t(1), dt, t(end));

    % Update max and min angular velocities if needed
    wMax = max(wMax, max(output{k}(:,5:7),[],'all'));
    wMin = min(wMin, min(output{k}(:,5:7),[],'all'));
end

% Add plotting buffer
wMax = wMax*buffer; 
wMin = wMin*buffer;

%% Plot results
for k = 1:length(w0)
    time = output{k}(:,1);
    X = output{k}(:,2:7);

    plotTitle = sprintf("EOM Evolution vs. Time for \\omega_0 = ^B[%.1f; %.1f; %.1f]", w0{k}(1), w0{k}(2), w0{k}(3));

    figure
    sgtitle(plotTitle)
    
    subplot(2,3,1)
    hold on
    plot(time, X(:,1));
    ylim([-1 1])
    xlabel("Time")
    ylabel("\sigma_1")
    
    subplot(2,3,2)
    hold on
    title("MRP Evolution")
    plot(time, X(:,2));
    ylim([-1 1])
    xlabel("Time")
    ylabel("\sigma_2")
    
    subplot(2,3,3)
    hold on
    plot(time, X(:,3));
    ylim([-1 1])
    xlabel("Time")
    ylabel("\sigma_3")
    
    subplot(2,3,4)
    hold on
    plot(time, X(:,4));
    ylim([wMin, wMax])
    xlabel("Time")
    ylabel("\omega_1")
    
    subplot(2,3,5)
    hold on
    title("Angular Velocity Evolution")
    plot(time, X(:,5));
    ylim([wMin, wMax])
    xlabel("Time")
    ylabel("\omega_2")
    
    subplot(2,3,6)
    hold on
    plot(time, X(:,6));
    ylim([wMin, wMax])
    xlabel("Time")
    ylabel("\omega_3")

end

