%% ASEN 5090 HW 4 Main Script
% By: Ian Faber, 10/03/2025

%% Housekeeping
clc; clear; close all;

addpath("..\utilities\");
addpath("..\HW3\HW3_Code");
addpath("..\HW3\HW3_Data\");

%% Setup
    % Choose a satellite
PRN = 14; % 2, 13, 16, 19, 20, 22 don't have valid C2L. 32 makes eph2pvt angry

    % User location
NIST_ECEF = [-1288398.567; -4721696.932; 4078625.350]; % m

    % Read RINEX and Ephemeris data
rinexData = rinexread("NIST00USA_R_20252230000_01D_30S_MO.rnx");
GPSData = rinexData.GPS;
PRNRnxIdx = GPSData.SatelliteID == PRN;
GPSData = GPSData(PRNRnxIdx, :);

ephData = read_clean_GPSbroadcast("brdc2230.25n", true);
PRNEphIdx = ephData(:,1) == PRN;
ephData = ephData(PRNEphIdx, :);

    % Planetary constants structure
pConst.c = 299792458; % m/s
pConst.wE = 7.2921151467e-5; % rad/s

%% Part 1. Compute expected range for each time from the RINEX file
    % Compute expected range
expRange = computeExpectedRange(GPSData, ephData, NIST_ECEF, pConst);

    % Extract C/A pseudorange measurement
C1C = GPSData.C1C;

    % Calculate pseudorange difference 0
dPR0 = C1C - expRange;

fprintf("dPR0 = C1C - expRange:\n")
fprintf("\tInitial dPR0 at %s = %.3f m\n", GPSData.Time(1), dPR0(1));
fprintf("\tFinal   dPR0 at %s = %.3f m\n\n", GPSData.Time(end), dPR0(end));

    % Plot difference
figure; tl = tiledlayout(1,3);
title(tl, "dPR0 = C1C - expRange")
nexttile([1 2]);
    hold on; grid on;
    title("dPR0 vs. time")
    plot(GPSData.Time, dPR0, 'b.')
    xlabel("Time [GPS]"); ylabel("dPR0 [m]"); yLim = ylim;
nexttile;
    hold on; grid on;
    title("dPR0 distribution")
    histogram(dPR0, 'FaceColor', 'b', 'Orientation', 'horizontal')
    xlabel("Bin Count"); ylim(yLim);

%% Part 2. Plot satellite clock correction in meters for each RINEX time
    % Extract time and calculate WN and TOW
t = GPSData.Time;
gpsEpoch = datetime(1980, 1, 6, 0, 0, 0);
tDiff = t - gpsEpoch;
WN = floor(seconds(tDiff)/(7*24*60*60));
TOW = mod(seconds(tDiff), 7*24*60*60);

    % Extract and plot clock correction and satellite position
[~, ephPos, ~, clkCorr, ~, ~] = eph2pvt2025(ephData, [WN, TOW], PRN);

figure;
hold on; grid on;
title("Clock Correction vs. Time")
plot(GPSData.Time, clkCorr, 'b.')
xlabel("Time [GPS]"); ylabel("Clock Correction [m]");

    % Calculate pseudorange difference 1
dPR1 = C1C - (expRange - clkCorr);

fprintf("dPR1 = C1C - (expRange - clkCorr):\n")
fprintf("\tInitial dPR1 at %s = %.3f m\n", GPSData.Time(1), dPR1(1));
fprintf("\tFinal   dPR1 at %s = %.3f m\n\n", GPSData.Time(end), dPR1(end));

    % Plot difference
figure; tl = tiledlayout(1,3);
title(tl, "dPR1 = C1C - (expRange - clkCorr)")
nexttile([1 2])
    hold on; grid on;
    title("dPR1 vs. time")
    plot(GPSData.Time, dPR1, 'b.')
    xlabel("Time [GPS]"); ylabel("dPR1 [m]"); yLim = ylim;
nexttile;
    hold on; grid on;
    title("dPR1 distribution")
    histogram(dPR1, 'FaceColor', 'b', 'Orientation', 'horizontal')
    xlabel("Bin Count"); ylim(yLim);

%% Part 3. Take relativity correction into account
    % Compute relativistic correction in meters
[~, ~, ~, ~, relCorr, ~, ~] = eph2pvtWithRel(ephData, [WN, TOW], PRN);

    % Calculate pseudorange difference 2
dPR2 = C1C - (expRange - clkCorr - relCorr);

fprintf("dPR2 = C1C - (expRange - clkCorr - relCorr):\n")
fprintf("\tInitial dPR2 at %s = %.3f m\n", GPSData.Time(1), dPR2(1));
fprintf("\tFinal   dPR2 at %s = %.3f m\n\n", GPSData.Time(end), dPR2(end));

figure;
hold on; grid on;
title("Relativistic Correction vs. Time")
plot(GPSData.Time, relCorr, 'b.')
xlabel("Time [GPS]"); ylabel("Relativistic Correction [m]");

    % Plot difference
figure; tl = tiledlayout(1,3);
title(tl, "dPR2 = C1C - (expRange - clkCorr - relCorr)")
nexttile([1 2])
    hold on; grid on;
    title("dPR2 vs. time")
    plot(GPSData.Time, dPR2, 'b.')
    xlabel("Time [GPS]"); ylabel("dPR2 [m]"); yLim = ylim;
nexttile;
    hold on; grid on;
    title("dPR2 distribution")
    histogram(dPR2, 'FaceColor', 'b', 'Orientation', 'horizontal')
    xlabel("Bin Count"); ylim(yLim);

%% Part 4. Compute simple tropospheric delay model and correct for it
    % Assumed zenith delay for NIST
zd = 2; % m

    % Calculate satellite elevation over this time interval
[~, satEl, ~] = compute_azelrange(NIST_ECEF, ephPos');

    % Calculate simple troposphere model
elevationCutoff = 5; % Ignore elevations below this angle
[tropoCorr, idxValid] = tropomodel(zd, deg2rad(satEl'), elevationCutoff);

    % Plot tropospheric correction
figure;
hold on; grid on;
title("Tropospheric Delay vs. time")
plot(GPSData.Time, tropoCorr, 'b.');
xlabel("Time [GPS]"); ylabel("Tropospheric Correction [m]")

    % Calculate pseudorange difference 3
dPR3 = C1C - (expRange - clkCorr - relCorr + tropoCorr);

%idxValid = ~isnan(dPR3);
validTime = GPSData.Time(idxValid);
validdPR3 = dPR3(idxValid);

fprintf("dPR3 = C1C - (expRange - clkCorr - relCorr + tropoCorr):\n")
fprintf("\tInitial dPR3 at %s = %.3f m, elevation cutoff = %.0f deg\n", validTime(1), validdPR3(1), elevationCutoff);
fprintf("\tFinal   dPR3 at %s = %.3f m, elevation cutoff = %.0f deg\n\n", validTime(end), validdPR3(end), elevationCutoff);

    % Plot difference
figure; tl = tiledlayout(1,3);
title(tl, "dPR3 = C1C - (expRange - clkCorr - relCorr + tropoCorr)")
nexttile([1 2])
    hold on; grid on;
    title("dPR3 vs. time")
    plot(GPSData.Time, dPR3, 'b.')
    xlabel("Time [GPS]"); ylabel(sprintf("dPR3 [m], \\theta_{el} < %.0f^o", elevationCutoff)); yLim = ylim;
nexttile;
    hold on; grid on;
    title("dPR3 distribution")
    histogram(dPR3, 'FaceColor', 'b', 'Orientation', 'horizontal')
    xlabel("Bin Count"); ylim(yLim);

%% Part 5. Compute ionospheric correction using an iono-free combination
    % Define carrier frequencies
f1 = 1575.42e6; % Hz
f2 = 1227.60e6; % Hz
f5 = 1176.45e6; % Hz

    % Extract C2L and C5Q signals
C2L = GPSData.C2L;
C5Q = GPSData.C5Q;

    % Compute ionospheric correction between L1 and L2 pseudoranges
[PRIF_12, iono_12] = ionocorr(C1C, f1, C2L, f2);

    % Plot ionospheric correction for each carrier frequency
figure;
hold on; grid on;
title("Ionospheric Correction for L1 and L2")
plot(GPSData.Time, iono_12(:,1), 'b.', GPSData.Time, iono_12(:,2), 'r.');
xlabel("Time [GPS]"); ylabel("Ionospheric Correction [m]");
legend("L1", "L2")

    % Calculate pseudorange difference 4
dPR4 = PRIF_12 - (expRange - clkCorr - relCorr + tropoCorr);

idxValid = ~isnan(dPR4);
validTime = GPSData.Time(idxValid);
validdPR4 = dPR4(idxValid);

fprintf("dPR4 = PRIF_12 - (expRange - clkCorr - relCorr + tropoCorr):\n")
fprintf("\tInitial dPR4 at %s = %.3f m, elevation cutoff = %.0f deg\n", validTime(1), validdPR4(1), elevationCutoff);
fprintf("\tFinal   dPR4 at %s = %.3f m, elevation cutoff = %.0f deg\n\n", validTime(end), validdPR4(end), elevationCutoff);

    % Plot difference
figure; tl = tiledlayout(1,3);
title(tl, "dPR4 = PRIF_{12} - (expRange - clkCorr - relCorr + tropoCorr)")
nexttile([1 2])
    hold on; grid on;
    title("dPR4 vs. time")
    plot(GPSData.Time, dPR4, 'b.')
    xlabel("Time [GPS]"); ylabel(sprintf("dPR4 [m], \\theta_{el} < %.0f^o", elevationCutoff)); yLim = ylim;
nexttile;
    hold on; grid on;
    title("dPR4 distribution")
    histogram(dPR4, 'FaceColor', 'b', 'Orientation', 'horizontal')
    xlabel("Bin Count"); ylim(yLim);

%% Part 6. Plot dPR1, dPR2, dPR3, and dPR4 on the same plot vs. time
figure;
grid on; hold on;
title("\rho Differences After Error Corrections")
plot(GPSData.Time, dPR1, 'b.');
plot(GPSData.Time, dPR2, 'r.');
plot(GPSData.Time, dPR3, 'k.');
plot(GPSData.Time, dPR4, 'm.');
xlabel("Time [GPS]"); ylabel("dPR [m]");
legend("dPR1 - satellite clock error", "dPR2 - relativity", "dPR3 - troposphere", "dPR4 - iono-free PR (using L1 and L2)", 'Location', 'best')

%% Part 7. Investigate pseudorange multipath
    % Extract signals of interest
C1C = GPSData.C1C;
L1C = GPSData.L1C;
S1C = GPSData.S1C;

C2W = GPSData.C2W;
L2W = GPSData.L2W;
S2W = GPSData.S2W;

C2L = GPSData.C2L;
L2L = GPSData.L2L;
S2L = GPSData.S2L;

C5Q = GPSData.C5Q;
L5Q = GPSData.L5Q;
S5Q = GPSData.S5Q;

    % Calculate multipath and CMC for C1C, C2W, C2L, and C5Q
[MP_C1C, CMC_C1C] = mpath(C1C, L1C, f1, L2L, f2);
[MP_C2W, CMC_C2W] = mpath(C2W, L2W, f2, L1C, f1);
[MP_C2L, CMC_C2L] = mpath(C2L, L2L, f2, L1C, f1);
[MP_C5Q, CMC_C5Q] = mpath(C5Q, L5Q, f5, L1C, f1);

    % Plot CMC and SNR vs. time
figure; tl = tiledlayout(2,1); ax = [];
title(tl, "CMC and SNR vs. Time")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    title("CMC vs. Time")
    plot(GPSData.Time, CMC_C1C, 'b-');
    plot(GPSData.Time, CMC_C2W, 'r-');
    plot(GPSData.Time, CMC_C2L, 'k-');
    plot(GPSData.Time, CMC_C5Q, 'm-');
    xlabel("Time [GPS]"); ylabel("CMC [m]");
    legend("C1C", "C2W", "C2L", "C5Q", 'Location', 'best');
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    title("SNR vs. Time")
    plot(GPSData.Time, S1C, 'b-');
    plot(GPSData.Time, S2W, 'r-');
    plot(GPSData.Time, S2L, 'k-');
    plot(GPSData.Time, S5Q, 'm-');
    xlabel("Time [GPS]"); ylabel("SNR [dB-Hz]");
    legend("S1C", "S2W", "S2L", "S5Q", 'Location', 'best');
linkaxes(ax, 'x');

    % Plot MP vs. time
figure;
hold on; grid on;
title("Multipath Observable vs. time")
plot(GPSData.Time, MP_C1C, 'b-');
plot(GPSData.Time, MP_C2W, 'r-');
plot(GPSData.Time, MP_C2L, 'k-');
plot(GPSData.Time, MP_C5Q, 'm-');
xlabel("Time [GPS]"); ylabel("Multipath Observable [m]")
legend("C1C w/ L2L", "C2W w/ L1C", "C2L w/ L1C", "C5Q w/ L1C", 'Location', 'best');

