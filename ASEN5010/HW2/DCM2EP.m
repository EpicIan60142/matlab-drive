%% ASEN 5010 DCM to Euler Parameters Script
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

%% Convert DCM to EP's

% CC 5/6 problem 4
theta = deg2rad([20, 10, -10]);

C = M1(theta(3))*M2(theta(2))*M3(theta(1)); % 3-2-1 rotation

% Sheppard's method
q0 = sqrt(0.25*(1 + trace(C)));
q1 = sqrt(0.25*(1 - trace(C) + 2*C(1,1)));
q2 = sqrt(0.25*(1 - trace(C) + 2*C(2,2)));
q3 = sqrt(0.25*(1 - trace(C) + 2*C(3,3)));

q = [q0, q1, q2, q3];

[max, idx] = max(q);

if idx == 1 % q0 was largest, can divide by it safely
    q1 = (C(2,3)-C(3,2))/(4*q0);
    q2 = (C(3,1)-C(1,3))/(4*q0);
    q3 = (C(1,2)-C(2,1))/(4*q0);
    q0 = (C(2,3)-C(3,2))/(4*q1); % Check sign on q0
elseif idx == 2 % q1 was largest, can divide by it safely
    q0 = (C(2,3)-C(3,2))/(4*q1);
    q2 = (C(1,2)+C(2,1))/(4*q1);
    q3 = (C(3,1)+C(1,3))/(4*q1);
    q1 = (C(2,3)-C(3,2))/(4*q0); % Check sign on q1
elseif idx == 3 % q2 was largest, can divide by it safely
    q0 = (C(3,1)-C(1,3))/(4*q2);
    q1 = (C(1,2)+C(2,1))/(4*q2);
    q3 = (C(2,3)+C(3,2))/(4*q2);
    q2 = (C(3,1)-C(1,3))/(4*q0); % Check sign on q2
else
    q0 = (C(1,2)-C(2,1))/(4*q3);
    q1 = (C(3,1)+C(1,3))/(4*q3);
    q2 = (C(2,3)+C(3,2))/(4*q3);
    q3 = (C(1,2)-C(2,1))/(4*q0); % Check sign on q3
end

EP = [q0, q1, q2, q3]'

