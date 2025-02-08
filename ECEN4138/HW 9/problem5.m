%% ECEN 4138 HW 9 Problem 6.7b
% - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
s = tf('s');

L = (s+2)/(s^2*(s+10)*(s^2+6*s+25));

%% Make Bode Plot
bode(L);
