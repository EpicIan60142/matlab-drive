%% ECEN 4138 HW 11 Problem 6.49
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
s = tf('s');

t = 0:0.001:100;

K = 5;

G = K/(s*(s/5 + 1)*(s/50 + 1));
C = (s+0.2)/(s+0.01);

L = C*G;

figure
margin(L)

figure
step(feedback(L,1),t)

figure
lsim(feedback(L,1),t,t)
