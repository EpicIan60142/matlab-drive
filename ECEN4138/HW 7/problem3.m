%% ECEN 4138 HW 7 Problem 5.11a
% - Ian Faber

%% Housekeeping
clc; clear; close all;

addpath('..\')

%% Define L(s)

s = tf('s');

L1 = (s^2+s+2)/(s*(s+5)*(s+6)*(s^2+2*s+1)); % L(s) for 5.11a

L = L1; % Choose root locus

%% Plot root locus
rlocus(L)
