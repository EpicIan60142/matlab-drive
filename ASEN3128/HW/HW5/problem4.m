%% Housekeeping
clc; clear; close all;

%% Constants
u0 = 502; % ft/s
K1 = -5.0;
K2 = -0.005;

%% Short Period Approx.
A_sp = [
            -1.019, 454.21;
            -0.005, -1.38
        ];

[V1, D1] = eig(A_sp)

%% With control
A_c = [
            -1.0189855, 446.91;
            -0.00500199, -0.38
        ];

[V2, D2] = eig(A_c)

%% Full state space model

A_full = [
            -0.02, 0.016, -0.65, -32.17;
            -0.13, -1.019, 454.21, 0;
            0, -0.005, -1.38, 0;
            0, 0, 1, 0
         ];

B = [-0.244; -1.46; -0.2; 0];

K = -[0 K2/u0 K1 0];

A_aug = A_full - B*K;

[V3, D3] = eig(A_aug)
