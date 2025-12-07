%% ASEN 5090 HW 2 Problem 4 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

addpath("HW2_Code\")
addpath("HW2_DATA\")

addpath("../utilities/")

%% Setup
NIST_ECEF = [-1288398.567, -4721696.932, 4078625.35]; % Given in ECEF
NIST_LLA = ecef2lla(NIST_ECEF);
Smead_LLA = [40.010245, -105.243979, 1601.07]; % Given in LLA
Smead_ECEF = lla2ecef(Smead_LLA);
EQUA_LLA = [0, NIST_LLA(2), 0];
EQUA_ECEF = lla2ecef(EQUA_LLA);

sp3Data = read_sp3("IGS0OPSFIN_20252230000_01D_15M_ORB.SP3");
satECEF = sp3Data(:,4:6)'*1000;
PRNs = sp3Data(:,3);

userPositions = [NIST_ECEF; Smead_ECEF; EQUA_ECEF]';

names = ["NIST", "Smead Backyard", "EQUA"];

%% Part a-c. Calculate azimuth and elevation angles for all satellites 
% throughout the day and make skyplots for each location

for k = 1:size(userPositions, 2)
        % Choose the correct user position
    userECEF = userPositions(:,k);
    
        % Calculate measurements
    [satAz, satEl, ~] = compute_azelrange(userECEF, satECEF);
    
    for kk = 1:length(satEl)
        if satEl(kk) < 0
            satAz(kk) = NaN;
            satEl(kk) = NaN;
        end
    end

        % Make sky plot
    fig = figure;
    title(sprintf("Skyplot for %s", names(k)))
    plotAzEl(satAz, satEl, PRNs, fig);
    % hold off;

end


