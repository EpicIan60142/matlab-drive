%% ASEN 5067 Final Project Controller Design Script
%   By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup plant

% Rotor inertia modeling
d = 0.004; % m
h = 0.011; % m

rho = 2710; % kg/m^3, assuming aluminum rotor

M = pi*(d/2)^2*h*rho; % kg

Jm = 0.5*M*(d/2)^2; % kg*m^2, rotor rotational inertia
JL = 1.7558e-6; % kg*m^2, Load inertia (3D printed arrow out of PLA, inertia from Autodesk Inventor)
J = Jm + JL;

% Motor specs
Ra = 10.0; % Ohm
Kt = 0.1226; % Nm/A
Ke = Kt; % V/rad/s

%% Design controller
Kp = 1.5;
Kd = 0.05;
Ki = 0.75;

s = tf('s');

G = (Kt/Ra)/(s*(J*s + (Kt*Ke)/Ra));
C = Kd*s + Kp + Ki/s;
H = 1/(s+3424.658); % Encoder pole, sends a sample every 292 us

Gmin = minreal(G)

% L = feedback(C*G,1)/Kp;
L = C*G;
T = feedback(L,1);

Lmin = minreal(L)

figure
hold on
margin(L)

figure
hold on
rlocus(L)

figure
hold on
step(T,0:0.001:10)


