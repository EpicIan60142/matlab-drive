%% ASEN 5090 HW 7 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

addpath("HW7_DATA\")
addpath("../utilities/")
addpath("../HW2/HW2_Code/")
addpath("../HW6/")

%% Setup
    % Read almanac data
almanacData = read_GPSyuma("YUMA236.alm", 2);

    % Read collected data
load("HW7data.mat");
s_R = signal;
clear signal

    % Approximate data collection location
userECEF = lla2ecef([40, -105.244, 1620])';

    % Time of collection in [WN, TOW] format
tCollect = datetime(2025,8,24,22,39,12,'TimeZone','UTC');
GPSEpoch = datetime(1980,1,6,0,0,0,'TimeZone','UTC');
dt = between(GPSEpoch, tCollect, {'weeks','time'});

WN = calweeks(dt);
TOW = seconds(time(dt));

    % Frequencies
Fs = 8e6; % Nominal sampling frequency
fIF = 20e3; % Nominal intermediate frequency

%% Part 1: Determine satellite visibility
    % Extract PRNs from the almanac file
PRNs = almanacData(:,1);

    % Populate positions
satECEF = [];
for k = 1:length(PRNs)
    [~, satPos, ~] = alm2pos(almanacData, [WN, TOW], PRNs(k));
    satECEF = [satECEF, satPos'];
end

    % Calculate azimuth, elevation, and range
[satAz, satEl, satRange] = compute_azelrange(userECEF,satECEF);

for kk = 1:length(satEl)
    if satEl(kk) < 0
        satAz(kk) = NaN;
        satEl(kk) = NaN;
    end
end

    % Make sky plot
fig = figure;
title(sprintf("Skyplot at data location, %s UTC", tCollect))
plotAzEl(satAz, satEl, PRNs, fig, true);

%% Part 2: Carrier wipeoff and code correlation for PRN 3
    % Assign PRN
PRN = 3;

    % Create time vector for 1 ms of data
tEnd = 1e-3;
dTs = 1/Fs;
tSample = 0:dTs:tEnd-dTs;

    % Create C/A code for desired PRN
        % Generate code sequence
fCA = 1.023e6; % C/A code frequency
G1 = ones(1,10);
G2 = ones(1,10);
CA = generatePRN(G1, G2, PRN, 1023);
    
        % Convert code to time domain
            % Populate time for C/A code itself
tCA = zeros(1, length(CA));
for k = 0:length(CA)
    tCA(k+1) = (1/fCA)*k;
end

            % Find code repeat period
tRepeat = length(CA)/fCA;

            % Map code to sample times
timeIdx = 1;
value = CA(timeIdx);
CA_PRN3 = zeros(1,length(tSample));
for k = 1:length(tSample)
    if tSample(k) >= tCA(timeIdx)
        value = CA(timeIdx);
        timeIdx = mod(timeIdx + 1, length(CA));
            % If the code repeats, increment the current time by the repeat
            % period
        if timeIdx == 0
            timeIdx = 1;
            tCA = tCA + tRepeat;
        end
    end
    CA_PRN3(k) = value;
end

first10 = CA(1:10)

CA_PRN3 = convPRNZeroOne2PosNeg(CA_PRN3);

    % Create doppler shift and code delay test values, then populate phases
fD = [0; 3.5e3];
shift = [9; 294];

theta = 2*pi*(fIF + fD).*tSample; % First row is fD = 0, second row is fD = 1e3 Hz

    % Test implementation
S_1 = 0;
for k = 1:length(tSample)
    S_1 = S_1 + s_R(k + shift(1))*CA_PRN3(k)*exp(-1j*theta(1,k));
end
fprintf("S_1 = %.4f + %.4fj\n\n", real(S_1), imag(S_1));

S_2 = 0;
for k = 1:length(tSample)
    S_2 = S_2 + s_R(k + shift(2))*CA_PRN3(k)*exp(-1j*theta(2,k));
end
fprintf("S_2 = %.4f + %.4fj\n\n", real(S_2), imag(S_2));

%% Part 3: Create a search grid in delay and doppler to find PRN 3
delay = 0:1023;
doppler = -20e3:500:20e3;
PRN = 3;

[DOP, TAU] = meshgrid(doppler, delay);

S = zeros(length(delay), length(doppler));

fprintf("Searching for PRN %.0f:\n", PRN);
for k = 1:length(delay)
    if ~mod(k-1,100) && k ~= 1
        fprintf("\tDelay = %.0f\n", delay(k)) % Print progress every 100 delays
    end
    for kk = 1:length(doppler)
        theta = 2*pi*(fIF + doppler(kk)).*tSample;
        for t = 1:length(tSample)
            S(k,kk) = S(k,kk) + s_R(t + delay(k))*CA_PRN3(t)*exp(-1j*theta(t));
        end
        S(k,kk) = norm(S(k,kk));
    end
end

[maxAmp, maxIdx] = max(S, [], 'all');
S = S/maxAmp;

fig = figure;
hold on; grid on;
title("Search Grid Results for PRN 3, 1 ms integration time")
res = surf(DOP, TAU, S, 'EdgeColor', 'none');
datatip(res, 'DataIndex', maxIdx);
xlabel("Doppler Shift [Hz]"); ylabel("Delay [Chips]"); zlabel("$\frac{|S(f,\tau)|}{|S_{max}|}$", 'interpreter', 'latex')
view([30 35]); colormap('cool')




