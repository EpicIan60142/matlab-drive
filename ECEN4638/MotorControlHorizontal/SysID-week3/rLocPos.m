%% ECEN 4638 Root Locus Script
%
% Ian Faber, Brennen Billig, Luke Hanley
%
%% Housekeeping
clc; clear; close all

%% Load Data
load("motorSystem_Data.mat");

%% Run rltool
modelsys = tf([4.9], [0.0886 1 0]);

rltool(modelsys)