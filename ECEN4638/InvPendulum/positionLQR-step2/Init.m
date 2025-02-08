%% ECEN 4638 Initialize Script
%
% Ian Faber, Brennen Billig, Luke Hanley
%
%% Housekeeping
clc; clear; close all;

%% Setup
Ts = 0.002; % Sampling time
Te = 20; % Experiment window
Tt = 5; % Transient Time

wc = 150; % Rate filter cutoff frequency
wr = 2*pi/Te; % Frequency resolution (smallest detectable frequency)
wn = pi/Ts; % Nyquist frequency (largest frequency to avoid aliasing)

N = 600; % Number of sinusoids to feed the model for SysID

wf = [1:N]*wr;
phaseRand = 2*pi*rand(1,N);

syms a b e d f

M = [1 0 0 0; 0 a 0 -b; 0 0 1 0; 0 -b 0 a];
Ahat = [0 1 0 0; 
        0 -e 0 0; 
        0 0 0 1; 
        0 0 d 0];
Bhat = [0; f; 0; 0];

A = M^-1*Ahat
B = M^-1*Bhat

% Set variables
Jeq = 0.0035842; % kgm^2
m = 0.125; % kg
r = 0.215; % m
etag = 0.9; % n.d.
Kg = 70; % n.d.
Jm = 3.87e-7; % kgm^2
L = 0.17675; % m
g = 9.81; % m/s^2
Beq = 0.004; % N/m/s
etam = 0.69; % n.d.
Kt = 0.00767; %Nm/A
Km = 0.00767; % rad/s/V
Rm = 2.6; % ohm

a = Jeq + m*r^2 + (etag*(Kg^2)*Jm);
b = m*L*r;
c = (4/3)*m*L^2;
d = m*g*L;
e = Beq + (etam*etag*Kt*Kg^2*Km)/Rm;
f = (etam*etag*Kt*Kg)/Rm;

M = [1 0 0 0; 0 a 0 -b; 0 0 1 0; 0 -b 0 c];
Ahat = [0 1 0 0; 
        0 -e 0 0; 
        0 0 0 1; 
        0 0 d 0];
Bhat = [0; f; 0; 0];

A = M^-1*Ahat
B = M^-1*Bhat

C = eye(4)
D = zeros(4,1)

modelSys = ss(A, B, eye(4), zeros(4,1))

%% Make LQR
un = 0.075; 
alphan = deg2rad(15);
thetan = pi/8; % rad
wthetan = 2000; % rad
walphan = 2000; % rad

xn = [thetan, alphan, wthetan, walphan];

Q = diag(1./xn.^2);
R = 1/un^2;

[Kc, P, clPoles] = lqr(modelSys, Q, R);
Kc


