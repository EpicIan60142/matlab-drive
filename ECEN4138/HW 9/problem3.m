%% ECEN 4138 HW 9 Problem 6.5d
% - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
s = tf('s');

L = (s^2+1)/(s*(s^2+4));

%% Make Bode Plot
bode(L);
