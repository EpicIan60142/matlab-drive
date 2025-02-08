%% ASEN 5014 HW 4
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Problem 2
M = [
          3 -1   4  0   7
          3  7  11 -9   8
          1 -3  -1  3   2
        -10  6 -11 -3 -23
    ]

r = rank(M);

[Q, R, P] = qr(M);

Q1 = Q(:,1:r)
Q2 = Q(:,r+1:end)

[Qhat, Rhat, Phat] = qr(M');

Qhat1 = Qhat(:,1:r)
Qhat2 = Qhat(:,r+1:end)

