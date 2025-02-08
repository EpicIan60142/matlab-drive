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

wr = 100; % Rate filter frequency

%% State space
A = [
        0   1;
        0   -w0
    ]

B = [
        0;
        K*w0
    ]

C = eye(2)
D = zeros(size(B))

Kp = 10 % Determined by not saturating voltage input

w_cl = sqrt(Kp*K*w0)
zeta = 0.6901 % 5% overshoot

Kd = (2*zeta*w_cl - w0)/(K*w0)

Kc = [Kp Kd];


