%% ECEN 4138 HW 11 Problem 6.58
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
s = tf('s');

z = 46.985; % rad/s

t = 0:0.001:100;

G = 30.287/(s*(s/20 + 1)*((s^2/100^2) + s/200 + 1));
C = (s + z)/(s + 3.16);

L = C*G;

figure
margin(L)

figure
step(feedback(L,1),t)

figure
lsim(feedback(L,1),t,t)
