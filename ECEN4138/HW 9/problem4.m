%% ECEN 4138 HW 9 Problem 6.6d
% - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
s = tf('s');

L = (s+3)/(s^2*(s+10));

%% Make Bode Plot
bode(L);