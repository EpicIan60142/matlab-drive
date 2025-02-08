%% ECEN 4138 HW 9 Problem 6.9
% - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
s = tf('s');

L = (s+10)/(s);

%% Make Bode Plot
figure
bode(L);

figure
step(L)
xlim([0, 3])