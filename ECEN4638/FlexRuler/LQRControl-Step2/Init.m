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

wc = 150; % Rate filter cutoff frequency
wr = 2*pi/Te; % Frequency resolution (smallest detectable frequency)
wn = pi/Ts; % Nyquist frequency (largest frequency to avoid aliasing)

N = 600; % Number of sinusoids to feed the model for SysID

wf = [1:N]*wr;
phaseRand = 2*pi*rand(1,N);

%% Set up state space
load("flexRulerSysID.mat", "T", "model1","model2")

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

%% Set up LQR
thetan = 1; % rad, step command
epsilonn = 1; % unitless, don't want too much strain
wn = 1000; % rad/s, want theta to respond quickly
epsilonDotn = 100; % rad/s, don't want strain to be too fast

un = 5; % V, don't saturate command input

xn = [thetan, epsilonn, wn, epsilonDotn];

Q = diag(1./xn.^2)
R = 1/un^2

[Kc, P, clPoles] = lqr(Ahat, Bhat, Q, R)

%% Make notch filter
Q = 3;
w1 = 95.504; % rad/s
w2 = 126.92; % rad/s

notch1 = (s^2 + w1^2)/(s^2 + (w1/Q)*s + w1^2)
notch2 = (s^2 + w2^2)/(s^2 + (w2/Q)*s + w2^2)

notch = notch1*notch2;

% figure
% bode(notch)

[num, den] = butter(10,w1,'s');
lpf = tf(num,den);
figure
margin(Kc*notch*trueSys)

