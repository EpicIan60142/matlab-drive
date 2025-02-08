%% ECEN 4638 Initialize Script
%
% Ian Faber, Brennen Billig, Luke Hanley
%
%% Housekeeping
clc; clear; close all;

%% Setup
Ts = 0.01; % Sampling time
T = 10; % Measurement window
Tt = 1; % Transient Time

wr = 2*pi/T; % Frequency resolution/smallest frequency
wn = pi/Ts; % Nyquist frequency/largest frequency to avoid aliasing

N = 200; % Number of sinusoids to feed the model
