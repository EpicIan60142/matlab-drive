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

wc = 150; % theta rate filter cutoff frequency
wc2 = 5; % ball rate filter cutoff frequency
wr = 2*pi/Te; % Frequency resolution (smallest detectable frequency)
wn = pi/Ts; % Nyquist frequency (largest frequency to avoid aliasing)

N = 200; % Number of sinusoids to feed the model for SysID

wf = [1:N]*wr;
phaseRand = 2*pi*rand(1,N);

%% Design inner loop PD

epsilon = 0.1;

KpInner = 8;
KdInner = 0.5/48;

%% Design outer loop PD

KpOuter = 0.9/9;
KdOuter = 0.3;
KiOuter = 0.05;

