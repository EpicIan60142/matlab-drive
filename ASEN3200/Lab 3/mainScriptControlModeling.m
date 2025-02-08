%% ASEN 3200 Lab A-3 Control Script
% Group 25 - Ian Faber, John Davis, Steven Liu, Ben Slama

%% Housekeeping
clc; clear; close all;

%% Constants

% System constants
K1 = 77.59/1000; % mNm -> Nm
K2 = 27.16/1000; % mNm -> Nm
I = 0.0068; % kgm^2

% Command pulse constants
period = 10;
offset = 0.72;
numPulses = 6;
stepCommand = 0.5;

% Transfer function constants
s = tf('s');
wn = sqrt(K1/I);
zeta = (K2/(2*I))*sqrt(I/K1);

% Poles
p1 = -zeta*wn + wn*sqrt(1-zeta^2)*1i;
p2 = -zeta*wn - wn*sqrt(1-zeta^2)*1i;

% Transfer functions
H1 = (wn^2)/(s^2 + 2*zeta*wn*s + wn^2);
% H2 = wn^2/((s-p1)*(s-p2));

%% Data Extraction
controlData = readmatrix("Data\task4Data");

time = (controlData(2:end,1) - controlData(2,1))/1000; % ms to sec
refPos = controlData(2:end,2); % rad
measPos = controlData(2:end,3)+0.01; % rad
current = controlData(2:end,4); % A
Kp = controlData(2:end,5); 
Kd = controlData(2:end,6);
Ki = controlData(2:end,7);

%% Analysis
t = 0:0.03:108;

pulse = 0;%((t > offset & t < period + offset) | (t > 2*period + offset & t < 3*period + offset) | (t > 4*period + offset & t < 5*period + offset) | (t > 6*period + offset & t < 7*period + offset) | (t > 8*period + offset & t < 9*period + offset) | (t > 10*period + offset & t < 11*period + offset))

for k = 0:numPulses-1
    pulse = pulse | (t > 2*k*period + offset & t < (2*k+1)*period + offset);
end

command = @(t) stepCommand*pulse;

[outputModel, timeModel] = lsim(H1, command(t), t);
% [outputTest, timeTest] = lsim(H2, command(t), t);

%% Plotting

figure
hold on
titleStr = sprintf("Step Response with K_p = %.2f mNm/rad, K_d = %.2f mNm/rad/s", mean(Kp), mean(Kd));
title(titleStr)
plot(time, measPos, 'b-')
plot(timeModel, outputModel, 'g-')
% plot(timeTest, outputTest, 'm-')
plot(time, refPos, 'r--')
xlabel("Time (sec)")
ylabel("Position (rad)")

legend("Step Response - Data", "Step Response - Model", "Step Command")

figure
hold on
titleStr = sprintf("Modeled Step Response with K_p = %.2f mNm/rad, K_d = %.2f mNm/rad/s", K1*1000, K2*1000);
title(titleStr)
plot(timeModel, outputModel, 'b-')
yline(1.05*stepCommand, 'k--')
yline(0.95*stepCommand, 'k--')
xlim([0, 10])
xlabel("Time (sec)")
ylabel("Position (rad)")

figure
hold on
titleStr = sprintf("Current Supplied to Reaction Wheel during Step Response");
title(titleStr)
yyaxis left
plot(time, current, 'b-');
ylabel("Current (A)")
yyaxis right
plot(time, refPos, 'r--')
ylabel("Commanded position (rad)")
ylim([-1 1])
xlabel("Time (sec)")

legend("Current Response", "Step Command")
