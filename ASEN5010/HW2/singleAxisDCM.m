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

%% Math

v = [3;7;2];

%% Problem CC3p2
theta = deg2rad([20,-10,120]);
C = M1(theta(3))*M2(theta(2))*M3(theta(1)) % 3-2-1 rotation

Phi = acos(0.5*(trace(C)-1))
e = (1/(2*sin(Phi)))*[C(2,3)-C(3,2); C(3,1)-C(1,3); C(1,2)-C(2,1)]



