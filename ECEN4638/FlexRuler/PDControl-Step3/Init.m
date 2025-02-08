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

wc = 20; % Rate filter cutoff frequency
wr = 2*pi/Te; % Frequency resolution (smallest detectable frequency)
wn = pi/Ts; % Nyquist frequency (largest frequency to avoid aliasing)

N = 600; % Number of sinusoids to feed the model for SysID

wf = [1:N]*wr;
phaseRand = 2*pi*rand(1,N);

%% Set up state space
load("flexRulerSysID.mat", "T", "model1","model2","model3")

s = tf('s');

oldSys = ss(T,'minimal')

A = oldSys.A;
B = oldSys.B;
C = oldSys.C;
D = oldSys.D;

Ahat = C*A*C^-1;
Bhat = C*B;
Chat = eye(4);
Dhat = zeros(4,1);

newSys = ss(Ahat, Bhat, Chat, Dhat) % Open Loop TF, states are [theta, epsilon, omega, and epsilonDot]

T = [model1; model2; s*model1; s*model2];
trueSys = ss(T,'minimal')
Atrue = trueSys.A
Btrue = trueSys.B
Ctrue = trueSys.C
Dtrue = trueSys.D

% bode(trueSys)

%% Set up PD
Kp = 15;
Kd = 1;

Kep = 3;

C = Kp*(1+(Kd/Kp)*s) + Kep;

%% Set up LPF
wfilt = 100; % rad/s
wfilt2 = 2;

lpf = 1/((s/wfilt) + 1);

lpf2 = 1/((s/wfilt2) + 1);






