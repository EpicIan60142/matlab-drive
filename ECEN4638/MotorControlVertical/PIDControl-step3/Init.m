%% ECEN 4638 Initialize Script
%
% Ian Faber, Brennen Billig, Luke Hanley
%
%% Housekeeping
clc; clear; close all;

%% Setup
load("Integrator-LeadData.mat")

Ts = 0.002; % Sampling time
Te = 15; % Experiment window
Tt = 5; % Transient Time

wc = 100; % Rate filter cutoff frequency
wr = 2*pi/Te; % Frequency resolution (smallest detectable frequency)
wn = pi/Ts; % Nyquist frequency (largest frequency to avoid aliasing)

N = 200; % Number of sinusoids to feed the model for SysID

K = 1.427;
w0 = 6.28319;
zeta = 1.0637;
Q = 1/(2*zeta);
s = tf('s');

H = K/((s^2/w0^2) + (1/(Q*w0))*s + 1);

% rltool(H)