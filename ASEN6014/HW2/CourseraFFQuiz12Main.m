%% ASEN 6014 Coursera Formation Flying Quiz 12
% By: Ian Faber, 09/26/2025

%% Housekeeping
clc; clear; close all;

addpath("../utilities/")

%% Setup
    % Gravitational parameter for Earth
mu = 398600.4415; % km^3/s^2

    % Provided orbital elements
oe.a = 7500; % km
theta = deg2rad(13);
oe.i = deg2rad(22);
q1 = 0.00707107;
q2 = 0.00707107;
oe.argPeri = atan2(q2, q1);
oe.e = sqrt(q1^2 + q2^2);
oe.f = theta - oe.argPeri;
oe.RAAN = deg2rad(70);

%% Create A matrix
A = calcCart2RelOrbA(oe, mu);

fprintf("\tCalculated A matrix:\n\n")

for k = 1:6
    fprintf("row%.0f = [%.6e,%.6e,%.6e,%.6e,%.6e,%.6e]\n", k, A(k,:));
end

%% Create A inverse matrix
Ainv = calcCart2RelOrbAinv(oe, mu);

fprintf("\n\tCalculated Ainv matrix:\n\n")

for k = 1:6
    fprintf("row%.0f = [%.6e,%.6e,%.6e,%.6e,%.6e,%.6e]\n", k, Ainv(k,:));
end

%% Check answer
check = A*Ainv - eye(6)

