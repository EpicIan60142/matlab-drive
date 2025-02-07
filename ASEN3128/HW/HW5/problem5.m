%% Housekeeping
clc; clear; close all;

%% A matrix
A = [
        -0.045, 0.036, 0, -32.2;
        -0.369, -2.02, 176, 0;
        0.0019, -0.0396, -2.948, 0;
        0, 0, 1, 0
    ];

[V, D] = eig(A)