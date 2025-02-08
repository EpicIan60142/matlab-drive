%% ECEN 4638 Initialize Script
%
% Ian Faber, Brennen Billig, Luke Hanley
%
%% Housekeeping
clc; clear; close all;

%% Setup
Ts = 0.002; % Sampling time
Te = 10; % Experiment window
Tt = 5; % Transient Time

wc = 200; % Rate filter cutoff frequency
wr = 2*pi/Te; % Frequency resolution (smallest detectable frequency)
wn = pi/Ts; % Nyquist frequency (largest frequency to avoid aliasing)

N = 200; % Number of sinusoids to feed the model for SysID

wf = [1:N]*wr;
phaseRand = 2*pi*rand(1,N);

%% Set up state space
load("flexJointSysID.mat", "T")

oldSys = ss(T,'minimal')

A = oldSys.A;
B = oldSys.B;
C = oldSys.C;
D = oldSys.D;

Ahat = C*A*C^-1;
Bhat = C*B;
Chat = eye(4);
Dhat = zeros(4,1);

newSys = ss(Ahat, Bhat, Chat, Dhat) % Open Loop TF, states are [theta1, theta2, w1, w2]

%% Set up LQR
theta1n = 1; % rad, step command
theta2n = 0.5; % rad, don't want theta2 to change
w1n = 1000; % rad/s, want theta1 to respond quickly
w2n = 100; % rad/s, don't want w2 to be too fast
en = 2; % Error, rad

un = 10; % V, don't saturate command input

xn = [theta1n, theta2n, w1n, w2n, en];

Q = diag(1./xn.^2)
R = 1/un^2

Chat2 = [1, 0, 0, 0];
Dhat2 = 0;

newSys = ss(Ahat, Bhat, Chat2, Dhat2);

[Kc, P, clPoles] = lqi(newSys, Q, R)



