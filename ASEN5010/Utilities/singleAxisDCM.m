%% ASEN 5010 Single-Axis DCM Script
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Definitions

M1 = @(theta) [
                1,          0,          0;
                0,  cos(theta), sin(theta);
                0,  -sin(theta), cos(theta)
              ];

M2 = @(theta) [
                cos(theta), 0,  -sin(theta);
                0,          1,            0;
                sin(theta), 0,  cos(theta)
              ];


M3 = @(theta) [
                cos(theta), sin(theta), 0;
                -sin(theta), cos(theta), 0;
                0,          0,          1
              ];


