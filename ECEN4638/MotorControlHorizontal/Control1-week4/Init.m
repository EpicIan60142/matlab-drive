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

wr = 2*pi/Te; % Frequency resolution (smallest detectable frequency)
wn = pi/Ts; % Nyquist frequency (largest frequency to avoid aliasing)

N = 200; % Number of sinusoids to feed the model for SysID

K = 4.4649;
w0 = 11.629;
s = tf('s');

H = K/(s*(s/w0 + 1));

% rltool(H)

%controller
C = 3.699;
z = w0;
p = 31.41;