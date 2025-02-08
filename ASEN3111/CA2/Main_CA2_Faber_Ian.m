%% ASEN 3111 - Computational Assignment 02 - Main
%   Driver program for CA 02 - Flow Over Thin Airfoils
%       Program that plots the flow over a thin airfoil using superposition
%       and the vortex panel method
%
%   Author: Ian Faber
%   Collaborators: Maggie Wussow, Ryan Caputo, Alex McCulley, Anna Sophia
%                  Rorrer Warren
%   Date Started: 2/23/2023
%   Date Finished: 2/28/2023

%% Housekeeping
clc; clear; close all;

tic;

%% Starter Command Window Management
fprintf("\t+--------------------------------+\n")
fprintf("\t\tCA 2 Results - Ian Faber\n\n\n")

%% Constants
c = 5; % m
alpha = 15; % deg
V_inf = 34; % m/s
p_inf = 101.3e3; % Pa = N/m^2
rho_inf = 1.225; % kg/m^3

bounds = [-2, 7; -3, 3]; % Min/max coordinates: [xmin, xmax; ymin, ymax]
nPoints = 100; % Number of points for grid (nPoints x nPoints)

N = 1000; % Number of vortices
NMax = 500; % Max number of vortices for error part
nLevels = [50, 100, 100]; % Streamlines, equipotential lines, pressure contours 

%% Bullet point 1
disp("Begin bullet point 1")

Plot_Airfoil_Flow(c, alpha, V_inf, p_inf, rho_inf, N, nLevels, bounds, nPoints);

drawnow % Plot the results for bullet point 1 immediately
toc; % Find time elapsed for bullet point 1;

%% Bullet point 2
disp("Begin bullet point 2")

Vortex_Error(c, alpha, V_inf, NMax, bounds, nPoints);

drawnow % Plot the results for bullet point 2 immediately
toc; % Find time elapsed up through bullet point 2

%% Bullet point 3
disp("Begin bullet point 3")

% Define study space for bullet point 3 (6 datapoints)
cStudy = 1:6;
alphaStudy = 0:3:15;
V_infStudy = 10:5:35;

% Adjust levels for visibility during the study
nLevelsStudy = [25,35,50];

Sensitivity_Analysis(cStudy, alphaStudy, V_infStudy, p_inf, rho_inf, N, nLevelsStudy, bounds, nPoints)

%% End Code
toc; % Find total elapsed time

%% Final Command Window Management
fprintf("\n\n\t+--------------------------------+\n\n")
disp(getRalphie)
