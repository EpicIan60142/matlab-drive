%% ASEN 6014 HW 1 Coursera Quiz 10 Question 7 Main Script
% By: Ian Faber, 08/28/2025

%% Housekeeping
clc; clear; close all;

addpath('../utilities')

%% Setup
cart.mu = 398600.4415; % km^3/s^2
cart.rVec = [-820.865; -1905.95; -7445.9]; % km
cart.vVec = [-6.75764; -1.85916; 0.930651]; % km/s

%% Convert cartesian to classical OE's
oe = convCart2ClassicOE(cart)

% Coursera answer for convenience
sma = [7500];
ecc = [0.05]; 
inc = [1.7802];
AN = [-2.8449];
AP = [2.6180];
f = [2.2690];
