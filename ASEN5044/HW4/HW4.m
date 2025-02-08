%% ASEN 5044 HW 4 Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Problem 4 - Simon 2.16
% Setup
min = -0.5;
max = 0.5;
nSamples = 10000;
nBins = 50;

% Generate random variables
x1 = min + (max-min)*rand(nSamples, 1);
x2 = min + (max-min)*rand(nSamples, 1);
x3 = min + (max-min)*rand(nSamples, 1);
x4 = min + (max-min)*rand(nSamples, 1);

% Plot histograms
A = (x1 + x2)/2;
B = (x1 + x2 + x3 + x4)/4;

figure(1)
hold on;
title("Histogram of (x_1 + x_2)/2")
hist(A, nBins)

figure(2)
hold on;
title("Histogram of (x_1 + x_2 + x_3 + x_4)/4")
hist(B, nBins)
