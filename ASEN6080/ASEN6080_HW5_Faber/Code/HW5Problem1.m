%% ASEN 6080 HW 5 Problem 1 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath("..\")
addpath(genpath("..\..\..\Utilities\"));

    % Extract parameters for Earth
earthConst = getPlanetConst();

    % Set up stations
stations = makeStations(earthConst);

    % Extract orbital elements
orbital = getOrbitConst();

    % ODE45 options
opt = odeset('AbsTol',1e-12,'RelTol',1e-12);

%% Problem 1a. Compare SRIF to LKF





