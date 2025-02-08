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
prefix = '..\Data\2022_08_25_013_Unit05_MEMS_';
% numFiles = 4;
% 
% for k = 1:numFiles
%     dataFiles{k} = dir([prefix,'*']);
% end

dataFiles = {[prefix, 'MAN'], [prefix, 'Fp2Cp5'], [prefix, 'Fp2C1'], [prefix, 'Fp4Cp5']};
clear prefix

%% Data extraction

% Trial 1: Frequency 0.2 Hz, Current 0.5 A
% Trial 2: Frequency 0.2 Hz, Current 1 A
% Trial 3: Frequency 0.4 Hz, Current 0.5 A

manualData = GyroData(dataFiles{1});
manualData.inputRate = manualData.inputRate * rpm2rad;

for k = 2:length(dataFiles)
    motorData{k-1} = GyroData(dataFiles{k});
    motorData{k-1}.inputRate = motorData{k-1}.inputRate * rpm2rad;
end

%% Data Analysis and plotting, part c.1

% Chosen trial: trial 1

% Compare gyro output to input rate
raw = figure(1);
hold on
grid on
title("Trial 1 Raw Gyro Output and Input Rate vs. Time")
xlabel("Time (sec)")
ylabel("Angular Rate (rad/s)")
plot(motorData{1}.time, motorData{1}.gyroOutput, '.');
plot(motorData{1}.time, -motorData{1}.gyroOutput, '.');
plot(motorData{1}.time, motorData{1}.inputRate, '.');
legend("Gyro Output", "Gyro Output (flipped)", "Input Rate", 'Location', 'best');

filename = "RawOutputVsTime";
saveas(raw, filename, imageFormat);

%% Data Analysis plotting, part c.2

% Chosen trial: trial 1

% Plot gyro output vs. input rate
cal1 = figure(2);
hold on
grid on
title("Trial 1 Gyro Output vs. Input Rate");
xlabel("Input Rate (rad/s)")
ylabel("Gyro Output (rad/s)")
plot(motorData{1}.inputRate, motorData{1}.gyroOutput, '.');

% Fit a line to figure 2's data and add it to the figure
[coef, approxCurve] = leastSquares(motorData{1}.inputRate, motorData{1}.gyroOutput, 1);

K(1) = coef(1);
bias(1) = coef(2);
label = sprintf("Line of best fit: \n gyroOutput = %.3f*inputRate + %.3f", K(1), bias(1));

t = linspace(min(motorData{1}.inputRate)-1, max(motorData{1}.inputRate)+1, 1000);

plot(t, approxCurve(t), 'k--', 'LineWidth', 1);
legend("Data", label, 'Location', 'best')

filename = "Trial1Calibration";
saveas(cal1, filename, imageFormat);

figure(3)
hold on
grid on
title("Trial 1 Adjusted Gyro Output and Input Rate vs. Time")
xlabel("Time (sec)")
ylabel("Angular Rate (rad/s)")

% Compute adjusted gyro output with calculated scale factor and bias
motorData{1}.adjustedGyroOutput = (1/K(1))*(motorData{1}.gyroOutput - bias(1));

plot(motorData{1}.time, motorData{1}.adjustedGyroOutput, '.');
plot(motorData{1}.time, motorData{1}.inputRate, '.');
legend("Adjusted Gyro Output", "Input Rate", 'Location', 'best')

%% Data Analysis and Plotting, part c.3

cal2 = figure(4);
hold on
grid on
title("Trial 2 Gyro Output vs. Input Rate");
xlabel("Input Rate (rad/s)")
ylabel("Gyro Output (rad/s)")
plot(motorData{2}.inputRate, motorData{2}.gyroOutput, '.');

[coef, approxCurve] = leastSquares(motorData{2}.inputRate, motorData{2}.gyroOutput, 1);
K(2) = coef(1);
bias(2) = coef(2);
label = sprintf("Line of best fit: \n gyroOutput = %.3f*inputRate + %.3f", K(2), bias(2));

t = linspace(min(motorData{2}.inputRate)-1, max(motorData{2}.inputRate)+1, 1000);

plot(t, approxCurve(t), 'k--', 'LineWidth', 1);
legend("Data", label, 'Location', 'best')

filename = "Trial2Calibration";
saveas(cal2, filename, imageFormat);


cal3 = figure(5);
hold on
grid on
title("Trial 3 Gyro Output vs. Input Rate");
xlabel("Input Rate (rad/s)")
ylabel("Gyro Output (rad/s)")
plot(motorData{3}.inputRate, motorData{3}.gyroOutput, '.');

[coef, approxCurve] = leastSquares(motorData{3}.inputRate, motorData{3}.gyroOutput, 1);
K(3) = coef(1);
bias(3) = coef(2);
label = sprintf("Line of best fit: \n gyroOutput = %.3f*inputRate + %.3f", K(3), bias(3));

t = linspace(min(motorData{3}.inputRate)-1, max(motorData{3}.inputRate)+1, 1000);

plot(t, approxCurve(t), 'k--', 'LineWidth', 1);
legend("Data", label, 'Location', 'best')

filename = "Trial3Calibration";
saveas(cal3, filename, imageFormat);

K = K';
bias = bias';
trial = (1:3)';

gyroTable = table(trial, K, bias);

KMean = mean(K);
Kstd = std(K);

biasMean = mean(bias);
biasSTD = std(bias);

%% Data Analysis and Plotting part c.4

% Chosen trials: 2 & 3

for k = 1:3
    % Compute adjusted gyro output with calculated scale factor and bias
    motorData{k}.adjustedGyroOutput = (1/K(k))*(motorData{k}.gyroOutput - bias(k));

    % Compute the error between the calibrated gyro output and the input
    % rate
    motorData{k}.rateError = motorData{k}.adjustedGyroOutput - motorData{k}.inputRate;
%     motorData{2}.time = motorData{2}.time(isfinite(motorData{2}.rateError));
%     motorData{2}.rateError = motorData{2}.rateError(isfinite(motorData{2}.rateError));
    
    % Calculate true and measured angular position vs. time
    motorData{k}.truePosition = cumtrapz(motorData{k}.time, motorData{k}.inputRate);
    motorData{k}.measuredPosition = cumtrapz(motorData{k}.time, motorData{k}.adjustedGyroOutput);

    % Calculate the error between the measured angular position and the
    % true angular position
    motorData{k}.positionError = motorData{k}.measuredPosition - motorData{k}.truePosition;
    
    % Clean up calculated values (remove Inf, if they exist) and rescale vectors
    valid = isfinite(motorData{k}.rateError) & isfinite(motorData{k}.positionError);
    
    motorData{k}.time = motorData{k}.time(valid);
    motorData{k}.inputRate = motorData{k}.inputRate(valid);
    motorData{k}.gyroOutput = motorData{k}.gyroOutput(valid);
    motorData{k}.adjustedGyroOutput = motorData{k}.adjustedGyroOutput(valid);
    motorData{k}.rateError = motorData{k}.rateError(valid);
    motorData{k}.positionError = motorData{k}.positionError(valid);
    motorData{k}.truePosition = motorData{k}.truePosition(valid);
    motorData{k}.measuredPosition = motorData{k}.measuredPosition(valid);
    
    motorData{k}.meanRateError = mean(motorData{k}.rateError);
    motorData{k}.stdRateError = std(motorData{k}.rateError);

    motorData{k}.meanPosError = mean(motorData{k}.positionError);
    motorData{k}.stdPosError = std(motorData{k}.positionError);
    
    meanRateError(k) = motorData{k}.meanRateError;
    stdRateError(k) = motorData{k}.stdRateError;

    meanPosError(k) = motorData{k}.meanPosError;
    stdPosError(k) = motorData{k}.stdPosError;

    rate = figure();
    hold on
    grid on
    titleText = sprintf("Trial %.0f Adjusted Gyro Output and Input Rate vs. Time", k);
    title(titleText);
    xlabel("Time (sec)")
    ylabel("Angular Rate (rad/s)")
    plot(motorData{k}.time, motorData{k}.adjustedGyroOutput, '.');
    plot(motorData{k}.time, motorData{k}.inputRate, '.');
    legend("Adjusted Gyro Output", "Input Rate", 'Location', 'best')
    
    filename = sprintf("Trial%.0fRatesVsTime", k);
    saveas(rate, filename, imageFormat)
    
    rateErr = figure();
    hold on
    grid on
    titleText = sprintf("Trial %.0f Angular Rate Measurement Error vs. Time", k);
    title(titleText)
    xlabel("Time (sec)")
    ylabel("Angular Rate Error (rad/s)")
    plot(motorData{k}.time, motorData{k}.rateError)

    filename = sprintf("Trial%.0fRateErrorVsTime", k);
    saveas(rateErr, filename, imageFormat)
    
    pos = figure();
    hold on
    grid on
    titleText = sprintf("Trial %.0f True and Measured Angular Positions vs. Time", k);
    title(titleText)
    xlabel("Time (sec)")
    ylabel("Angular Position (rad)")
    plot(motorData{k}.time, motorData{k}.truePosition);
    plot(motorData{k}.time, motorData{k}.measuredPosition);
    legend("True Position", "Measured Position")

    filename = sprintf("Trial%.0fPositionsVsTime", k);
    saveas(pos, filename, imageFormat)

    posErr = figure();
    hold on
    grid on
    titleText = sprintf("Trial %.0f Angular Position Error vs. Time", k);
    title(titleText)
    xlabel("Time (sec)")
    ylabel("Angular Position Error (rad)")
    plot(motorData{k}.time, motorData{k}.positionError);

    filename = sprintf("Trial%.0fPositionErrorVsTime", k);
    saveas(posErr, filename, imageFormat)

end

close all;

TotalMeanRateError = mean(meanRateError);
TotalStdRateError = mean(stdRateError);

TotalMeanPosError = mean(meanPosError);
TotalStdPosError = mean(stdPosError);

%% Manual experiment
manual = figure();
hold on
grid on
title("Manual Input rate vs. Gyro Output")
xlabel("Time")
ylabel("Angular rate (rad/s)")
plot(manualData.time, manualData.gyroOutput, '.')