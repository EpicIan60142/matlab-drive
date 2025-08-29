%% ASEN 6014 HW 1 Coursera Quiz 10 Question 6 Main Script
% By: Ian Faber, 08/28/2025

%% Housekeeping
clc; clear; close all;

addpath('../utilities')

%% Setup
oe.mu = 398600.4415; % km^3/s^2
oe.a = 8000; % km
oe.e = 0.1;
oe.i = deg2rad(30); % rad
oe.RAAN = deg2rad(145); % rad
oe.argPeri = deg2rad(120); % rad

M0 = deg2rad(10); % rad

n = sqrt(oe.mu/(oe.a^3)); % rad/s

%% Propagate mean anomaly forward by 1 hour and convert to cartesian
t_end = 3600; % sec

M = M0 + n*(t_end - 0);
E = convM2E(M, oe.e);
oe.f = convE2f(E,oe.e);

cart = convClassicOE2Cart(oe);

rVec = cart.rVec
vVec = cart.vVec

% Coursera answer for convenience
rVec_coursera = [-1264.6, 8013.8, -3371.2];
vVec_coursera = [-6.0396, -0.2044, 2.0967];


