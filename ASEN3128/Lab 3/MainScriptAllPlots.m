%% ASEN 3128 Lab 3 Main Script
% Section 011 - Ian Faber, Felix Evrard, Blake Hogen, Gabriel Agostine

%% Housekeeping
clc; clear; close all;

%% Constants
m = 0.068; % kg
g = 9.81; % m/s^2
d = 0.06; % m
Km = 0.0024; % N*m/N
Ix = 5.8*10^-5; % kgm^2
Iy = 7.2*10^-5; % kgm^2
Iz = 1*10^-4; % kgm^2

I = [Ix, 0, 0; 0, Iy, 0; 0, 0, Iz];

nu = 10^-3; % N/(m/s)^2
mu = 2*10^-6; % M*m/(rad/s)^2

%% ODE45 Problem 2

tspan = [0 10];

X0 = zeros(12, 1);
X0(1:3) = [0; 0; -5]; % x y z
X0(4:6) = deg2rad([0; 0; 0]); % phi theta psi
X0(7:9) = [0; 0; 0]; % u v w
X0(10:12) = [0.1; 0.1; 0.1]; % p q r

options = odeset('Events', @detectGround);

Fc = [0; 0; -m*g];
Gc = zeros(3,1);

[time, state] = ode45(@(t, var)AircraftEOM(t, var, g, m, nu, mu, Fc, Gc), tspan, X0, options);

aircraft_state_array = [
                            state(:,1)';
                            state(:,2)';
                            state(:,3)';
                            state(:,4)';
                            state(:,5)';
                            state(:,6)';
                            state(:,7)';
                            state(:,8)';
                            state(:,9)';
                            state(:,10)';
                            state(:,11)';
                            state(:,12)'
                       ];

control_input_array = [
                            ones(1,length(time))*Fc(3);
                            ones(1,length(time))*Gc(1);
                            ones(1,length(time))*Gc(2);
                            ones(1,length(time))*Gc(3)
                      ];

motorVector1 = [];

for k = 1:length(time)
    motorForces = ComputeMotorForces([0;0;control_input_array(1,k)], control_input_array(2:4, k), d, Km);
    motorVector1 = [motorVector1, motorForces];
end


fig = 1:6;
col = 'b-';

PlotAircraftSim(time, aircraft_state_array, control_input_array, fig, col);

%% ODE45 problem 3

% tspan = [0 10];
% 
% X0 = zeros(12, 1);
% X0(1:3) = [0; 0; -5]; % x y z
% X0(4:6) = deg2rad([0; 0; 0]); % phi theta psi
% X0(7:9) = [0; 0; 0]; % u v w
% X0(10:12) = [0.0; 0.0; 0.1]; % p q r
% 
% deltaFc = [0; 0; 0];
% deltaGc = zeros(3, 1);
% 
% [time, state] = ode45(@(t,var)QuadrotorEOM_Linearized(t, var, g, m, I, deltaFc, deltaGc), tspan, X0);
% 
% aircraft_state_array = [
%                             state(:,1)';
%                             state(:,2)';
%                             state(:,3)';
%                             state(:,4)';
%                             state(:,5)';
%                             state(:,6)';
%                             state(:,7)';
%                             state(:,8)';
%                             state(:,9)';
%                             state(:,10)';
%                             state(:,11)';
%                             state(:,12)'
%                        ];
% 
% control_input_array = [
%                             ones(1,length(time))*Fc(3);
%                             ones(1,length(time))*Gc(1);
%                             ones(1,length(time))*Gc(2);
%                             ones(1,length(time))*Gc(3)
%                       ];
% 
% fig = 1:6;
% col = 'r-';
% 
% PlotAircraftSim(time, aircraft_state_array, control_input_array, fig, col)

%% ODE45 Problem 6

tspan = [0 10];

X0 = zeros(12, 1);
X0(1:3) = [0; 0; -5]; % x y z
X0(4:6) = deg2rad([0; 0; 0]); % phi theta psi
X0(7:9) = [0; 0; 0]; % u v w
X0(10:12) = [0.1; 0.1; 0.1]; % p q r

[time2, state] = ode45(@(t,var)QuadrotorEOMwithControl(t, var, g, m, nu, mu), tspan, X0);

aircraft_state_array = [
                            state(:,1)';
                            state(:,2)';
                            state(:,3)';
                            state(:,4)';
                            state(:,5)';
                            state(:,6)';
                            state(:,7)';
                            state(:,8)';
                            state(:,9)';
                            state(:,10)';
                            state(:,11)';
                            state(:,12)'
                       ];

control_input_array = [];

for k = 1:length(time2)
    [Fc, Gc] = RotationDerivativeFeedback(aircraft_state_array(:,k), m, g);
    control_input_array = [control_input_array, [Fc(3); Gc(1); Gc(2); Gc(3)]];
end

motorVector2 = [];

for k = 1:length(time2)
    motorForces = ComputeMotorForces([0;0;control_input_array(1,k)], control_input_array(2:4, k), d, Km);
    motorVector2 = [motorVector2, motorForces];
end

% control_input_array = [
%                             Fc(3);
%                             Gc(1);
%                             Gc(2);
%                             Gc(3);
%                       ];

fig = 1:6;
col = 'r-';

PlotAircraftSim(time2, aircraft_state_array, control_input_array, fig, col)

% % Plot motor forces
figure
sgtitle("Drone Motor Forces")
subplot(4,1,1)
hold on;
title("f1 vs. time")
plot(time, motorVector1(1,:), 'b'); hold on;
plot(time2, motorVector2(1,:), 'r'); hold on;
ylabel("f1 (N)")
xlabel("time (sec)")

subplot(4,1,2)
hold on;
title("f2 vs. time")
plot(time, motorVector1(1,:), col); hold on;
plot(time2, motorVector2(2,:), col); hold on;
ylabel("f2 (N)")
xlabel("time (sec)")

subplot(4,1,3)
hold on;
title("f3 vs. time")
plot(time, motorVector1(1,:), col); hold on;
plot(time2, motorVector2(3,:), col); hold on;
ylabel("f3 (N)")
xlabel("time (sec)")

subplot(4,1,4)
hold on;
title("f4 vs. time")
plot(time, motorVector1(1,:), col); hold on;
plot(time2, motorVector2(4,:), col); hold on;
ylabel("f4 (N)")
xlabel("time (sec)")

