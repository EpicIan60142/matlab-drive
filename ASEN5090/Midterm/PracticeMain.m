%% ASEN 5090 Practice Midterm Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

addpath("utilities\")

%% Setup
RE = 6380; % km
c = 3e5; % km/s
mu = 398600; % km^3/s^2

L1 = 1575.42e6; % Hz

%% Problem 2
    % Setup
R_T_LLA = [40; 30; 100]; % deg, deg, m

format long;

    % Part a
R_T_ECEF = lla2ecef(R_T_LLA')'

    % Part b
C_ECEF2ENU = ECEF2ENU(R_T_LLA(1), R_T_LLA(2))

    % Part c and e
R_T_ENU = C_ECEF2ENU*R_T_ECEF;
R_A_ENU = R_T_ENU + [-900; 0; -100];
R_A_ECEF = C_ECEF2ENU'*R_A_ENU
[Az, El, Range] = compute_azelrange(R_A_ECEF, R_T_ECEF)

    % Part d
LOS = R_T_ENU - R_A_ENU

    % Part e



format default;