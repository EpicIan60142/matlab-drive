%% ECEN 4638 Parameter ID Script
%
% Ian Faber, Luke Hanley, Brennen Billig
%
%% Housekeeping
clc; clear; close all;

%%  Variable setup
K = 4.9;
w0 = 11.282;
tau = 1/w0;
w1 = 100;

%% Open simulink model
open("SRV02_demo.slx")


