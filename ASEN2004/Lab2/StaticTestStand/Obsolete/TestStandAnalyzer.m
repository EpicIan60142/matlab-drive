%% ASEN 2004 Lab 2 Static Test Stand Analyzer
% Section 301 - The Diabetics of Space Z
% Tanner Brummond, Joshua Camp, Ian Faber, Joshua Goodrich, Justin
% McGregor, Andrew Vo
% Code by: Ian Faber
% Last Modified: 3/14/22, 12:25 AM

%% Housekeeping and Setup
clc; clear; close all;

sampFreq = 1652; % Sampled at 1.652 kHz
maxTime = 5; % Maximum allowable time for data start and stop calcs
relError = 0.05; % Relative error between thrust averages for start and stop calcs
g0 = 9.81; % Freefall acceleration on Earth
mWater = 1; % Mass of water propellant (1000 g)

%% Data Extraction
numFiles = 16;

for i = 1:numFiles
    files{i} = dir(['Static Test Stand Data\Fixed Mass\LA_Test_FixedMass_Trial',num2str(i)]);
    files{i}.path = ['Static Test Stand Data\Fixed Mass\LA_Test_FixedMass_Trial',num2str(i)];
end

%% Data Processing
for i = 1:numFiles
    files{i}.data = load(files{i}.path);
    files{i}.thrust = files{i}.data(:,3);
    files{i}.peakThrust = max(files{i}.thrust);
    files{i}.time = (linspace(0,length(files{i}.thrust)/sampFreq,length(files{i}.thrust)))'; % Compute time based on sampling frequency
    
    % Calculate a 5-value moving average for detecting trend changes
    thrustAvg = (files{i}.thrust(1:end-4) + files{i}.thrust(2:end-3) + files{i}.thrust(3:end-2) + files{i}.thrust(4:end-1) + files{i}.thrust(5:end))/5;
    
    % Check for abrupt changes in the data and if time is less than 5 seconds
    conditionStart = ischange(thrustAvg) & files{i}.time(1:end-4) < maxTime;
    
    % Check to see if data changes by less than relError*100% and if time is less than 5 seconds
    conditionStop = thrustAvg(1:end-1)./thrustAvg(2:end) - 1 > relError & files{i}.time(1:end-5) < maxTime;
    
    start(i) = find(conditionStart, 1, 'first');
    stop(i) = find(conditionStop, 1, 'last');
    
    files{i}.thrustTime = files{i}.time(stop(i))- files{i}.time(start(i));
    
    files{i}.waterFlow = mWater/files{i}.thrustTime; % Calculate mass flow rate of the water
    files{i}.impulse = trapz(files{i}.time(start(i):stop(i)), files{i}.thrust(start(i):stop(i)) - files{i}.waterFlow.*files{i}.time(start(i):stop(i)));
    files{i}.isp = files{i}.impulse/(mWater*g0);
    
    isp(i,1) = files{i}.isp(end);
    SEMIspPart(i,1) = std(isp)/sqrt(i);
    peakThrust(i,1) = files{i}.peakThrust;
    thrustTime(i,1) = files{i}.thrustTime;
    
end

avgISP = mean(isp)
stdISP = std(isp)

avgPeakThrust = mean(peakThrust)
stdPeakThrust = std(peakThrust)

avgThrustTime = mean(thrustTime)
stdThrustTime = std(thrustTime)

SEMIspFull = stdISP/sqrt(numFiles)

%% Plotting
for i = 1:numFiles
    figure
    subplot(1,2,1)
    hold on;
    title("Raw Thrust Data")
    plot(files{i}.time, files{i}.thrust);
    xline(files{i}.time(start(i)));
    xline(files{i}.time(stop(i)));
    xlabel("Time (sec)")
    ylabel("Thrust (N)")
    hold off;

    subplot(1,2,2)
    hold on;
    title("Truncated Thrust Data")
    plot(files{i}.time(start(i):stop(i)), files{i}.thrust(start(i):stop(i)));
    xline(files{i}.time(start(i)));
    xline(files{i}.time(stop(i)));
    xlabel("Time (sec)")
    ylabel("Thrust (N)")
    hold off;
end

