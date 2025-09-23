%% ASEN 5090 HW 3 Problem 3 Main Script
% By: Ian Faber, 09/21/2025

%% Housekeeping
clc; clear; close all;

addpath("..\utilities\")
addpath("HW3_Code\");
addpath("HW3_Data\");
addpath("..\HW2\HW2_Code\");
addpath("..\HW2\HW2_DATA\");

%% Setup
    % Choose a satellite
PRN = 5;

    % NIST ECEF coords
NIST_ECEF = [-1288398.567; -4721696.932; 4078625.35]; % Given in ECEF

    % Speed of light
c = 299792458; % m/s

    % Angular velocity of earth
wE = 7.2921151467e-5; % rad/s

    % Convergence criterion
epsilon = 1e-3;

%% Part a. Compute range and elevation of satellite from NIST based on broadcast epemeris data
    % Extract ephemeris data and compute positions
ephData = read_clean_GPSbroadcast('brdc2230.25n', true);

    % Extract only the satellite of interest
PRNIdx = ephData(:,1) == PRN;

ephData = ephData(PRNIdx, :);

    % Extract times to compute positions at
rinexData = rinexread("NIST00USA_R_20252230000_01D_30S_MO.rnx");
GPSData = rinexData.GPS;
PRNIdx = GPSData.SatelliteID == PRN;
GPSData = GPSData(PRNIdx, :);
t = GPSData.Time;

    % Calculate week number and time of week
gpsEpoch = datetime(1980, 1, 6, 0, 0, 0);
tDiff = t - gpsEpoch;
WN = floor(seconds(tDiff)/(7*24*60*60));
TOW = mod(seconds(tDiff), 7*24*60*60);

    % Calculate GPS satellite position using ephemeris data
[~, ephPos, ephVel, ~, ~, ~] = eph2pvt2025(ephData, [WN, TOW], PRN);

    % Calculate range and elevation of satellite from NIST
ephPos = ephPos';
[~, satEl, satRange] = compute_azelrange(NIST_ECEF, ephPos);
satEl = satEl';
satRange = satRange';

    % Plot vs. time
fig = figure; tl = tiledlayout(2,1); ax = [];
title(tl, sprintf("Range and Elevation vs. time for PRN %.0f from NIST", PRN))
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, satRange, 'b.', 'DisplayName', 'Computed Range')
    xlabel("Time"); ylabel("Range [m]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, satEl, 'k.')
    xlabel("Time"); ylabel("Elevation [deg]")
linkaxes(ax, 'x')

%% Part b. Compute the expected range accounting for signal travel time and frame rotation
    % Loop through all measurement times
expRange = [];
for k = 1:length(satRange)
        % Step 1
    wn_r = WN(k);
    t_r = TOW(k);
    pos_r = ephPos(:, k);

        % Step 2
    R_geo = norm(pos_r - NIST_ECEF);

        % Step 3
    R_old = R_geo;
    t_t = t_r - R_old/c;

        % Step 4
    [~, pos_t, ephVel, ~, ~, ~] = eph2pvt2025(ephData, [wn_r, t_t], PRN);

        % Step 5
    phi = wE*(t_r - t_t);
    pos_r = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1]*pos_t';

        % Step 6
    R_new = norm(pos_r - NIST_ECEF);

    while(abs((R_new/R_old)-1) > epsilon)
        R_old = R_new;

            % Step 3
        t_t = t_r - R_old/c;

            % Step 4
        [~, pos_t, ephVel, ~, ~, ~] = eph2pvt2025(ephData, [wn_r, t_t], PRN);
    
            % Step 5
        phi = wE*(t_r - t_t);
        pos_r = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1]*pos_t';
    
            % Step 6
        R_new = norm(pos_r - NIST_ECEF);
    end

    expRange = [expRange; R_new];
end

    % Plot
fig;
nt = nexttile(1);
    hold on;
    plot(t, expRange, 'r.', 'MarkerSize', 4, 'DisplayName', 'Expected Range')
legend(nt, 'Location', 'bestoutside')    

figure; tl = tiledlayout(1,3);
title(tl, sprintf("Difference of Computed and Expected range for PRN %.0f", PRN))
nexttile([1 2]);
    hold on; grid on;
    title("Computed Range - Expected Range")
    plot(satRange - expRange, 'b.');
    xlabel("Time"); ylabel("\Delta\rho [m]")
nexttile;
    hold on;
    title("Range Difference Histogram")
    histogram(satRange - expRange, 'BinWidth', 10, 'FaceColor', 'b', 'Orientation', 'horizontal')
    xlabel("Bin Count")

%% Part c. Plot the C1C pseudoranges and Expected range vs. time and their difference

figure;
hold on; grid on;
title(sprintf("C1C and Expected Range for PRN %.0f vs. time", PRN))
plot(t, GPSData.C1C, 'b.')
plot(t, expRange, 'r.', 'MarkerSize', 4)
xlabel("Time"); ylabel("Range [m]")
legend("C1C Pseudorange", "Expected Range", 'Location', 'bestoutside');

figure; tl = tiledlayout(1,3);
title(tl, sprintf("Difference between C1C and Expected range for PRN %.0f", PRN))
nexttile([1 2]);
    hold on; grid on;
    title("Range difference vs. time")
    plot(t, GPSData.C1C - expRange, 'b.')
    xlabel("Time"); ylabel("\Delta\rho [m]")
nexttile;
    hold on;
    title("Range Difference Histogram")
    histogram(GPSData.C1C - expRange, 'FaceColor', 'b', 'Orientation', 'horizontal')
    xlabel("Bin Count")



