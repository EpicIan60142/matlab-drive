%% ASEN 6014 Formation Flying Coursera Quiz 3 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

addpath("..\utilities\")

%% Setup
format long

    % Constants
const.mu = 398600; % km^3/s^2

    % Chief spacecraft
r_c_N = [-4893.268; 3864.478; 3169.646]; % km
v_c_N = [-3.91386; -6.257673; 1.59797]; % km/s

    % Problem 1
r_d_N = [-4892.98; 3863.073; 3170.619]; % km
v_d_N = [-3.913302; -6.258661; 1.598199]; % km/s

    % Problem 2
rho_H = [-0.537; 1.221; 1.106]; % km
rhoPrime_H = [0.000486; 0.001158; 0.0005590]; % km/s

    % Problem 3
x = 10; y = 500; r = 7000; % km
xDot = 0.1; yDot = -0.1; rDot = 0.05; % km/s

    % Problem 4
dr = 10; s = 500; r = r; % km
drDot = 0.1; sDot = -0.1; rDot = rDot; % km/s

%% Problem 1: Convert inertial deputy state to Hill frame
x_c_N = [r_c_N; v_c_N];
x_d_N = [r_d_N; v_d_N];

x_d_H = convDeputyN2H(x_c_N, x_d_N, const)

%% Problem 2: Convert Hill frame deputy to Inertial
x_d_H = [rho_H; rhoPrime_H];

x_d_N = convDeputyH2N(x_c_N, x_d_H, const)

%% Problem 3: Convert rectilinear coordinates to curvilinear coordinates
x_d_RL = [x; y; xDot; yDot];
r_c = [r; rDot];

x_d_CL = convDeputyRL2CL(x_d_RL, r_c)

%% Problem 4: Convert curvilinear coordinates to rectilinar coordinates
x_d_CL = [dr; s; drDot; sDot];

x_d_RL = convDeputyCL2RL(x_d_CL, r_c)

format default
