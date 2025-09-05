%% ASEN 5090 HW 1 Problem 3 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

addpath("../utilities/")

%% Setup
chips = 1023;
PRNs = [19, 25, 5];
CAs = [];

%% Generate PRN C/A codes
for k = 1:length(PRNs)
    % Reset registers and extract PRN
    G1 = ones(1,10);
    G2 = ones(1,10);
    PRN = PRNs(k);

    fprintf("\nGenerating C/A code for PRN %.0f", PRN)

    % Generate 1023-bit C/A code
    [CA, G1_end, G2_end] = generatePRN(G1, G2, PRN, 1023);
    
    CAs = [CAs; CA];

    % Report first 16 chips
    first16 = binaryVectorToHex(CA(1:16));
    fprintf("\nFirst 16 bits for PRN %.0f are 0x%s\n", PRN, first16)

    % Report last 16 chips
    last16 = binaryVectorToHex(CA(end-15:end));
    fprintf("Last 16 bits for PRN %.0f are 0x%s\n", PRN, last16)
    
    % Plot first and last 16 chips
    fig = figure; tl = tiledlayout(2,1);
    title(tl, sprintf("First and last 16 chips of PRN %.0f", PRN))
    nexttile;
        hold on; grid on; xticks(1:16)
        title(sprintf("First 16 chips: 0x%s", first16))
        stairs(0:15, CA(1:16)); xlim([0, 16])
        xlabel("Epoch [n.d.]"); ylabel("Chip [n.d.]")
    nexttile;
        hold on; grid on; xticks((length(CA)-15):length(CA))
        title(sprintf("Last 16 chips: 0x%s", last16))
        stairs((length(CA)-15):length(CA), CA(end-15:end));
        xlabel("Epoch [n.d.]"); ylabel("Chip [n.d.]")

end



