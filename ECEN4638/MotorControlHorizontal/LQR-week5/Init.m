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

% N = 200; % Number of sinusoids to feed the model for SysID

K = 4.4649;
w0 = 11.629;

wc = 100; % Rate Derivative cutoff frequency

A = [
        0   1;
        0   -w0
    ];

B = [
        0;
        K*w0
    ];

C = eye(2);
D = zeros(size(B));

sys = ss(A,B,C,D);

x1n = 1; % rad
x2n = 1000; % rad/s
un = 10; % V

xn = [x1n, x2n];

Q = diag(1./xn.^2)
R = 1/un^2

[Kc, P, clPoles] = lqr(A, B, Q, R)

