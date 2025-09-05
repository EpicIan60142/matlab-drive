%% ASEN 5090 HW 1 Problem 4 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

addpath("../utilities/")

%% Setup
nShift = 0:1022; % Vector of chip shifts to plot
G1 = ones(1,10); % Initial G1 register state
G2 = ones(1,10); % Initial G2 register state
chips = 1023; % Size of C/A codes to generate

%% Part a. Plot Autocorrelation for PRN 19
    % Define PRN and C/A code
PRN_auto = 19;
CA_auto = generatePRN(G1, G2, PRN_auto, chips);
CA_auto = convPRNZeroOne2PosNeg(CA_auto);
    
    % Calculate auto-correlation
Rn_auto = zeros(size(nShift));
for k = 1:length(nShift)
    Rn_auto(k) = calcPRNAutoCorr(CA_auto, nShift(k));
end

    % Plot auto-correlation
figure;
hold on; grid on;
title(sprintf("Autocorrelation of PRN %.0f", PRN_auto))
plot(nShift, Rn_auto)
xlabel("Chip shift n"); ylabel(sprintf("R$^{(%.0f)}$(n)", PRN_auto), 'Interpreter','latex');

%% Part b. Plot Cyclic Cross-Correlation between PRN 19 and 200-chip delayed PRN 19
    % Setup 
PRN_cross = 19;
nDelay = 200;

    % Generate base and delayed PRNs
CA_cross = generatePRN(G1, G2, PRN_cross, chips);
CA_cross = convPRNZeroOne2PosNeg(CA_cross);
CA_crossDelayed = shiftCA(CA_cross, nDelay);

    % Calculate cross-correlation
Rn_cross = zeros(size(nShift));
for k = 1:length(nShift)
    Rn_cross(k) = calcPRNCrossCorr(CA_cross, CA_crossDelayed, nShift(k));
end

    % Plot cross-correlation
figure;
hold on; grid on;
title(sprintf("Cyclic cross-correlation of PRN %.0f with PRN %.0f delayed by %.0f chips", PRN_cross, PRN_cross, nDelay))
plot(nShift, Rn_cross)
xlabel("Chip shift n"); ylabel(sprintf("R$^{(%.0f, %.0f)}$(n)", PRN_cross, PRN_cross), 'Interpreter','latex');

%% Part c/d. Plot Cyclic Cross-Correlation between PRN 19 and PRN 25, PRN 5
    % Setup
PRN_cross_base = 19;
PRN_cross_test = [25, 5];

    % Generate base PRN code
CA_k = generatePRN(G1, G2, PRN_cross_base, chips);
CA_k = convPRNZeroOne2PosNeg(CA_k);

    % Loop over secondary PRNs
for k = 1:length(PRN_cross_test)
        % Generate secondary PRN code
    CA_l = generatePRN(G1, G2, PRN_cross_test(k), chips);
    CA_l = convPRNZeroOne2PosNeg(CA_l);

        % Calculate cross correlation
    Rn_cross_test = zeros(length(PRN_cross_test), length(nShift));
    for kk = 1:length(nShift)
        Rn_cross_test(k, kk) = calcPRNCrossCorr(CA_k, CA_l, nShift(kk));
    end

        % Plot cross correlation
    figure;
    hold on; grid on;
    title(sprintf("Cyclic cross-correlation of PRN %.0f with PRN %.0f", PRN_cross_base, PRN_cross_test(k)))
    plot(nShift, Rn_cross_test(k, :))
    xlabel("Chip shift n"); ylabel(sprintf("R$^{(%.0f, %.0f)}$(n)", PRN_cross_base, PRN_cross_test(k)), 'Interpreter','latex');

end

%% Part e. Sum together delayed codes of PRN 19, PRN 25, and PRN 5, then correlate with undelayed PRN 19
    % Setup
PRNs = [19, 25, 5];
delays = [350, 905, 75];

    % Generate clean PRN 19
CA_19 = generatePRN(G1, G2, 19, chips);
CA_19 = convPRNZeroOne2PosNeg(CA_19);

    % Loop over constituent PRNs
CA_sum = zeros(1, chips);
CA_i = zeros(length(PRNs), chips);
for k = 1:length(PRNs)
        % Calculate PRN
    CA_i(k,:) = generatePRN(G1, G2, PRNs(k), chips);
    CA_i(k,:) = convPRNZeroOne2PosNeg(CA_i(k,:));

        % Add to C/A sum
    CA_sum = CA_sum + CA_i(k,:);
end

    % Calculate cross-correlation between clean PRN 19 and summed signal
Rn_sum = zeros(size(nShift));
for k = 1:length(nShift)
    Rn_sum(k) = calcPRNCrossCorr(CA_19, CA_sum, nShift(k));
end

    % Plot cross-correlation
figure;
hold on; grid on;
title(sprintf("Cyclic cross-correlation of PRN %.0f with sum of delayed PRNs %.0f, %.0f, %.0f", 19, PRNs))
plot(nShift, Rn_sum)
xlabel("Chip shift n"); ylabel(sprintf("R$^{(%.0f, %.0f+%.0f+%.0f)}$(n)", 19, PRNs), 'Interpreter','latex');

%% Part f. Plot simulated noise and PRNs 19, 25, and 5
    % Setup
std = 4;
noise = std*randn(1, chips);
lims = [-3*std, 3*std];

    % Plot all signals together
fig = figure; tl = tiledlayout(4,1); ax = [];
title(tl, sprintf("Delayed PRNs %.0f, %.0f, %.0f, and noise", PRNs));
nt = nexttile; ax(end+1) = nt;
    hold on; grid on;
    plot(CA_i(1,:));
    xlabel("Chip"); ylabel("x_1"); ylim(lims);
nt = nexttile; ax(end+1) = nt;
    hold on; grid on;
    plot(CA_i(2,:));
    xlabel("Chip"); ylabel("x_2"); ylim(lims);
nt = nexttile; ax(end+1) = nt;
    hold on; grid on;
    plot(CA_i(3,:));
    xlabel("Chip"); ylabel("x_3"); ylim(lims);
nt = nexttile; ax(end+1) = nt;
    hold on; grid on;
    plot(noise);
    xlabel("Element"); ylabel("noise"); ylim(lims);
linkaxes(ax, 'x')

%% Part g. Sum together x1 - x3 and noise, then correlate with fresh PRN 19
    % Sum together the constituent signals
CA_noisy = sum(CA_i) + noise;

    % Calculate cross-correlation with clean PRN 19
Rn_noisy = zeros(size(nShift));
for k = 1:length(nShift)
    Rn_noisy(k) = calcPRNCrossCorr(CA_19, CA_noisy, nShift(k));
end

    % Plot cross-correlation
figure;
hold on; grid on;
title(sprintf("Cyclic cross-correlation of PRN %.0f with noisy signal", 19))
plot(nShift, Rn_noisy)
xlabel("Chip shift n"); ylabel(sprintf("R$^{(%.0f, noise)}$(n)", 19), 'Interpreter','latex');

