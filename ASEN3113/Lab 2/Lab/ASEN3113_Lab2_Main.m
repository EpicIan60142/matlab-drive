%% ASEN 3113 Lab 2 Code
%   - Alyxis Ellington, Ian Faber, Joe Machi, Alex McCulley, Austin Marx

%% Housekeeping
clc; clear; close all;

%% Constants
% Aluminum, Brass, Steel
rho = [2810, 8500, 8000]; % kg/m^3
c = [960, 380, 500]; % J/kgK
k = [130, 115, 16.2]; % W/mK
alpha = k./(rho.*c); % m^2/s

d = 1*0.0254; % in -> m

%% Data Extraction
files = ["Data\Brass_26V_245mA", "Data\Brass_29V_273mA","Data\Steel_21V_192mA"];
fileNum = 2;

input = extractData(files(fileNum), rho, c, k, alpha, d);

