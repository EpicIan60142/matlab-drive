%% ASEN 3200 Lab A-1 Gyro Script
% Section 013 Group 28: Ian Faber, Ian Mccomas, Ben Slama
% Descripton: Script that processes raw gyro data and measures the bias and
% adjusted scale factor for the MEMS gyro, calibrates the output, and
% calculates angular position and angular position error from the
% calibrated gyro output

%% Housekeeping
clc; clear; close all;

%% Constants and conversions
rpm2rad = pi/30;
imageFormat = 'jpg';

%% Data file setup
prefix = '..\Data\2022_08_25_013_Unit05_RWHEEL_';
% numFiles = 4;
% 
% for k = 1:numFiles
%     dataFiles{k} = dir([prefix,'*']);
% end

dataFiles = {[prefix, 'T4t5'], [prefix, 'T5t5'], [prefix, 'T10t5'], [prefix, 'T15t5'], [prefix, 'T20t5']};
clear prefix

%% Data extraction

% Trial 1: Torque 4 mNm, Time 5 sec
% Trial 2: Torque 5 mNm, Time 5 sec
% Trial 3: Torque 10 mNm, Time 5 sec
% Trial 4: Torque 15 mNm, Time 5 sec
% Trial 5: Torque 20 mNm, Time 5 sec

for k = 1:length(dataFiles)
    torqueData{k} = RWData(dataFiles{k});
    torqueData{k}.angV = torqueData{k}.angV * rpm2rad;
end

%% Data analysis and Plotting, part c

for k = 1:5
    start = find(torqueData{k}.angV > 1.75, 1, 'first');
    [~, stop] = max(torqueData{k}.angV);

    [coef, approxCurve] = leastSquares(torqueData{k}.time(start:stop), torqueData{k}.angV(start:stop), 1);
    alpha(k) = coef(1);
    t = linspace(torqueData{k}.time(start) - 0.5, torqueData{k}.time(stop) + 0.5, 1000);
    label = sprintf("Line of best fit: \\omega(t) = %.3f*t + %.3f \n \\alpha_{avg}: %.3f rad/s^2", coef(1), coef(2), coef(1));

    ang = figure();
    hold on
    grid on
    titleText = sprintf("Trial %.0f Angular Velocity vs. Time", k);
    title(titleText);
    xlabel("Time (sec)")
    ylabel("Angular velocity (rad/s)")
    plot(torqueData{k}.time, torqueData{k}.angV, '.');
    plot(t, approxCurve(t), 'k--');
    xline(torqueData{k}.time(start), 'g--')
    xline(torqueData{k}.time(stop), 'r--')
    legend("Angular Velocity", label, "Acceleration start", "Acceleration stop", 'Location', 'best')
    
    filename = sprintf("Trial%.0fAngularVelocityVsTime", k);
    saveas(ang, filename, imageFormat);

    % Torque = I*alpha, I = Torque/alpha
    I(k) = torqueData{k}.cmdTorque(end)/alpha(k);
end

close all;

meanI = mean(I)
stdI = std(I)

%% Reaction Wheel Capacity
aeroTorque = 10^-4; % Nm
limit = 4000*rpm2rad; % 4000 rpm to rad/s

t = meanI*limit/aeroTorque

capacity = aeroTorque*t






