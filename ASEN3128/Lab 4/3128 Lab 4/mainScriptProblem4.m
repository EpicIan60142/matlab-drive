%% ASEN 3128 Lab 4 Main Script
%   Section 011 - Ian Faber, Sam Mohnacs, Blake Wilson, Kyle Elligott

%% Housekeeping
clc; clear; close all;

%% Constants
g = 9.81; % m/s^2
m = 0.068; % kg
d = 0.06; % m
Km = 0.0024; % N*m/N
Ix = 5.8*10^-5; % kg*m^2
Iy = 7.2*10^-5; % kg*m^2
Iz = 1.0*10^-4; % kg*m^2
nu = 10^-3; % N/(m/2)^2
mu = 2*10^-6; % N*m/(rad/s)^2

%% ODE45 Problem 4

tspan = [0 10];

X0 = zeros(12, 1);
X0(1:3) = [0; 0; -5]; % x y z
X0(4:6) = deg2rad([5; 0; 0]); % phi theta psi
X0(7:9) = [0; 0; 0]; % u v w
X0(10:12) = [0.0; 0.0; 0.0]; % p q r

options = odeset('Events', @detectGround);

[time, state] = ode45(@(t, var)QuadrotorEOMwithControl(t, var, g, m, nu, mu), tspan, X0, options);

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

for k = 1:length(time)
    [Fc, Gc] = InnerLoopFeedback(aircraft_state_array(:,k));
    control_input_array = [control_input_array, [Fc(3); Gc(1); Gc(2); Gc(3)]];
end

motorVector = [];

for k = 1:length(time)
    motorForces = ComputeMotorForces([0;0;control_input_array(1,k)], control_input_array(2:4, k), d, Km);
    motorVector = [motorVector, motorForces];
end


fig = 1:6;
col = 'b-';

PlotAircraftSim(time, aircraft_state_array, control_input_array, fig, col);


