%% ASEN 5090 HW 2 Problem 2 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

addpath("HW2_Code\")
addpath("HW2_DATA\")

%% Setup
NIST_ECEF = [-1288398.567, -4721696.932, 4078625.35]; % Given in ECEF
Smead_LLA = [40.010245, -105.243979, 1601.07]; % Given in LLA

%% Part a. Create table of ECEF positions and lat/long/height values for 
%%% NIST, Smead backyard, and EQUA

    % Define location names
names = ["NIST", "Smead Backyard", "EQUA"];

    % Convert NIST ECEF to LLA
NIST_LLA = ecef2lla(NIST_ECEF, 'WGS84');

    % Convert Smead LLA to ECEF
Smead_ECEF = lla2ecef(Smead_LLA, 'WGS84');

    % Calculate EQUA ECEF given latitude of 0 deg, longitude same as
    % NIST, and height of roughly 
EQUA_LLA = [0, NIST_LLA(2), 0];
EQUA_ECEF = lla2ecef(EQUA_LLA);

    % Report results in a table
variables = ["Location", "X_ECEF [m]", "Y_ECEF [m]", "Z_ECEF [m]", "Latitude [deg N]", "Longitude [deg E]", "Height [m]"];
data = [
            names(1), round(NIST_ECEF, 0), round(NIST_LLA, 6);
            names(2), round(Smead_ECEF, 0), round(Smead_LLA, 6);
            names(3), round(EQUA_ECEF, 0), round(EQUA_LLA, 6)
       ];

result = table(data(:,1), data(:,2), data(:,3), data(:,4), data(:,5), data(:,6), data(:,7), 'VariableNames', variables)

%% Part b. Write a function to calculate the DCM from ECEF to ENU
    % Define unit vectors for NIST and EQUA
r_NIST = NIST_ECEF';
r_NIST = r_NIST/norm(r_NIST);
r_EQUA = EQUA_ECEF';
r_EQUA = r_EQUA/norm(r_EQUA);

    % Calculate C_ECEF2ENU for each location
C_ECEF2ENU_NIST = ECEF2ENU(NIST_LLA(1), NIST_LLA(2));
C_ECEF2ENU_EQUA = ECEF2ENU(EQUA_LLA(1), EQUA_LLA(2));

    % Apply the DCMs and interpret the result
r_NIST_ENU = C_ECEF2ENU_NIST*r_NIST
r_EQUA_ENU = C_ECEF2ENU_EQUA*r_EQUA



