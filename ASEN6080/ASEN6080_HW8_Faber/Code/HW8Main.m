%% ASEN 6080 HW 8 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
    % Path logistics
addpath("..\")
addpath(genpath("..\..\..\Utilities\"));

    % Extract Earth params
pConst = getPlanetConst();

    % Extract spacecraft params
scConst = getSCConst();

    % Extract orbit params
orbital = getOrbitConst();

    % Number of Monte Carlo runs
N = 1000;


