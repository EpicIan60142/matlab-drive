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

%% Problem CC7p1
theta = deg2rad([20,10,-10]);
DCM321_1 = M1(theta(3))*M2(theta(2))*M3(theta(1)) % 3-2-1 rotation

t1_1 = atan2(DCM321_1(3,1),-DCM321_1(3,2))
t2_1 = acos(DCM321_1(3,3))
t3_1 = atan2(DCM321_1(1,3),DCM321_1(2,3))

theta = [t1_1,t2_1,t3_1];
DCM313_1 = M3(theta(3))*M1(theta(2))*M3(theta(1)) % 3-1-3 rotation

%% Problem CC8p1
theta1 = deg2rad([10, 20, 30]); % 3-2-1
theta2 = deg2rad([40.6387, 35.5311, -36.0535]); % 3-1-3

theta = theta1;
DCM321 = M1(theta(3))*M2(theta(2))*M3(theta(1)) % 3-2-1 rotation

t1 = atan2d(DCM321(3,1),-DCM321(3,2))
t2 = acosd(DCM321(3,3))
t3 = atan2d(DCM321(1,3),DCM321(2,3))

theta = deg2rad([t1,t2,t3]);
DCM313 = M3(theta(3))*M1(theta(2))*M3(theta(1)) % 3-1-3 rotation

% theta = theta2;
% DCM313 = M3(theta(3))*M1(theta(2))*M3(theta(1)) % 3-1-3 rotation

%% Problem CC8p2
theta3 = deg2rad([10,20,30]); % 3-2-1
theta4 = deg2rad([-5,5,5]); % 3-2-1
theta5 = deg2rad([13.2251, 16.3676, 23.6181]); % 3-2-1

% In-class example
% theta3 = deg2rad([30,-45,60]); % 3-2-1
% theta4 = deg2rad([10,25,-15]); % 3-2-1
% theta5 = deg2rad([-0.933242, -72.3373, 79.9636]);



theta = theta3;
BN = M1(theta(3))*M2(theta(2))*M3(theta(1)); % 3-2-1
theta = theta4;
RN = M1(theta(3))*M2(theta(2))*M3(theta(1)); % 3-2-1

BR = BN*RN'

yaw = atan2d(BR(1,2),BR(1,1))
pitch = -asind(BR(1,3))
roll = atan2d(BR(2,3),BR(3,3))

theta = deg2rad([yaw, pitch, roll]);
DCMRel = M1(theta(3))*M2(theta(2))*M3(theta(1)) % 3-2-1

% theta = theta5;
% DCMRel = M1(theta(3))*M2(theta(2))*M3(theta(1)) % 3-2-1


