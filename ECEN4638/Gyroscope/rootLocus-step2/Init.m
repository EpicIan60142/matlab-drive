%% ECEN 4638 Initialize Script
%
% Ian Faber, Brennen Billig, Luke Hanley
%
%% Housekeeping
clc; clear; close all;

%% Setup
Ts = 0.002; % Sampling time
Te = 15; % Experiment window
Tt = 5; % Transient Time

wc = 150; % Rate filter cutoff frequency
wr = 2*pi/Te; % Frequency resolution (smallest detectable frequency)
wn = pi/Ts; % Nyquist frequency (largest frequency to avoid aliasing)

N = 200; % Number of sinusoids to feed the model for SysID

wf = [1:N]*wr;
phaseRand = 2*pi*rand(1,N);

modelsys = tf([db2mag(6.023)], [1/3.351 1 0]);
% rltool(modelsys)

%lead controller
s = tf("s");
controlla = (9.727*(s + 3.35))/(s + 11.43);

% Gyro gain
Gg = 1.7341*1.025;
deadzone = 0.02; % phi deadzone, offset from the start

%PD controller
kp = 10; % maxes u out at 10
kd = 1.2; % achieves 5% overshoot



