clc; clear; close all;

% xdot1 = 3x1 + u1 + u2
% xdot2 = 2x2 + u2 + d

A = [3 0; 0 2];
B = [1 1 0; 0 1 1];
C = eye(2);
D = zeros(size(B));

sys = ss(A, B, C, D)

tFinal = 10;
t = 0:0.001:tFinal;

u = [0*ones(size(t)); 1*ones(size(t)); 1*ones(size(t))];
x0 = [0; 0];

lsim(sys, u, t)


