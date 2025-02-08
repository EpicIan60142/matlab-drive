%% ASEN 5067 Final Project doodling
% - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup plant
drag = 0; % Drag enabled?

P = 2.5; % W
T = 0.5; % kg*cm
V = 12; % V

d = 0.004; % m
h = 0.011; % m

rho = 7800; % kg/m^3

I = P/V;
Ra = V/I;

T = 9.807*T/100; % kg*cm -> Nm
Kt = T/I;

M = pi*(d/2)^2*h*rho; % kg

Jm = 0.5*M*(d/2)^2; % kg*m^2, rotor rotational inertia
JL = 0.1; % Load inertia
J = Jm + JL;

% Ra = 57.6; % Ohms => 2.5W/12V = 0.2083333 A => 12/0.2083333 = 57.6
% Kt = 2.4; % Nm/A => 0.005 kg*m 
Ke = Kt; % V/rad/s
% Jm = 0.00111; % kg*m^2, assumed iron as shaft metal. Cylinder w/ 4mm diameter and 11 mm height, rho = 7800 kg/m^3
b = 0.001*drag; % Nm/rad/s

Kp = 1;
Kd = 5;
Ki = 10;

s = tf('s');

G = (Kt/Ra)/(s*(J*s + b + (Kt*Ke)/Ra));
C = Kd*s + Kp + Ki/s;

minreal(G)

% L = feedback(C*G,1)/Kp;
L = C*G/(1+C*G);

minreal(L)

figure
hold on
rlocus(L)
