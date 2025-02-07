%% Project 1: Calorimetry to determine an unknown sample
% By: Ian Faber (108577813) and Justin Travis (109468486)

% Created: 9/25/2021
% Last modified: 10/20/2021

% Purpose: Identify sample B by calculating the specific heat and
% associated uncertainty in context of 5 given materials.

% Inputs: Calorimetry temperature profile data, calorimeter specific heat,
% sample and calorimeter masses and uncertainties, and systematic
% thermocouple uncertainty.

% Outputs: Least squares regression slopes and intercepts, extrapolated 
% temperatures, sample specific heat, uncertainties for each variable, 
% overall specific heat uncertainty, range of specific heat, percent errors 
% from given materials.

%% Set up and preliminary plotting
clear; clc; close all;

%Read in data
data = readmatrix("Sample_B.txt");
data(:,1) = [];

%Start at the 3rd entry to avoid NaN values at the beginning
time = data(3:end,1);
waterTemp = data(3:end,2);
sampleTemp1 = data(3:end,3);
sampleTemp2 = data(3:end,4);

%Average two thermocouple readings, since errors will be practically
%identical since the same type of thermocouple is used for both
meanTemp = (sampleTemp1 + sampleTemp2)/2;

%Initialize extraneous variables
Cc = 0.895;
Ms = 39.3;
sigma_Ms = 0.001;
Mc = 510.0;
sigma_Mc = 0.05;
sigma_thermo = 2;

%Plot given data
figure

%Water temp readings
subplot(2,2,1);
hold on;
grid on;
waterPlot = plot(time, waterTemp, 'b-');
T1Temp = plot(time,mean(waterTemp)*ones(length(time)), 'r-');
ylim([92 95]);
title("Water temperature vs. time");
xlabel("time (sec)");
ylabel("Temperature (^oC)");
legend(T1Temp, "T1");

%Temperature profile readings
subplot(2,2,2);
hold on;
grid on;
plot(time, sampleTemp1, 'r-');
plot(time, sampleTemp2, 'g-');
plot(time, meanTemp, 'k-');
title("Calorimeter temperature vs. time");
xlabel("time (sec)");
ylabel("Temperature (^oC)");
legend("First sample reading", "Second sample reading", "Average sample reading",'Location','best')

%Estimate important time stamps from data by looking at the temperature curve and find the vector indices where they occur
timeAdded = 217.273; %Time sample is added
timeNoCurveLow = 276.583; %Time when the curve of the data starts noticeably
timeNoCurveHigh = 503.841; %Time when the curve of the data ends noticeably
sampleStart = find(time == timeAdded); %Index when sample is added
sampleNoCurveLow = find(time == timeNoCurveLow); %Index where curved part of the profile begins
sampleNoCurveHigh = find(time == timeNoCurveHigh); %Index where the curved part of the profile ends

%Split up time vector into useful sections
firstSection = time(1:sampleStart);
middleSection = time(sampleStart:sampleNoCurveLow);
lastSection = time(sampleNoCurveHigh:end);
curveSection = time(sampleNoCurveLow:sampleNoCurveHigh);

%Used to separate the curved section from the others, purely for
%visualization
offsetX = 0;
offsetY = 0;

%Plot different sections together
subplot(2,2,4)
hold on;
grid on;
firstPlot = plot(firstSection, meanTemp(1:sampleStart), 'k--');
middlePlot = plot(middleSection, meanTemp(sampleStart:sampleNoCurveLow),'k--');
lastPlot = plot(lastSection, meanTemp(sampleNoCurveHigh:end),'k--');
curvePlot = plot(curveSection+offsetX, meanTemp(sampleNoCurveLow:sampleNoCurveHigh)+offsetY,'k--');

%% Least Squares regressions and plotting

%Least squares on first section
d = meanTemp(1:sampleStart);
Acol1 = firstSection;
Acol2 = ones(length(firstSection),1);
A = [Acol1 Acol2];

x_hat_1 = linReg(A, d)
LeastSquares1 = @(x) x_hat_1(1)*x + x_hat_1(2); %Anonymous function handle for the first regression

%Plot regression line for first section
firstReg = plot(firstSection, LeastSquares1(firstSection));

%Least squares on second section
e = meanTemp(sampleStart:sampleNoCurveLow);
Bcol1 = middleSection;
Bcol2 = ones(length(middleSection),1);
B = [Bcol1 Bcol2];

x_hat_2 = linReg(B, e)
LeastSquares2 = @(x)x_hat_2(1)*x + x_hat_2(2); %Anonymous function handle for the second regression

%Least squares on third section
f = meanTemp(sampleNoCurveHigh:end);
Ccol1 = lastSection;
Ccol2 = ones(length(lastSection),1);
C = [Ccol1 Ccol2];

x_hat_3 = linReg(C, f)
LeastSquares3 = @(x)x_hat_3(1)*x + x_hat_3(2); %Anonymous function handle for the third regression

%Find the time and index where the second and third regressions intersect
intersectT = (x_hat_3(2)-x_hat_2(2))/(x_hat_2(1)-x_hat_3(1));

intersect = find((time < (intersectT)), 1, 'last' ); %Find max index

%Plot second and third regressions
secondReg = plot(time(sampleStart:intersect), LeastSquares2(time(sampleStart:intersect)));
lastReg = plot(time, LeastSquares3(time));

%% Temperature extrapolation and plotting

%Find T0, TH, and Tmid
T0 = LeastSquares1(time(sampleStart))
T1 = mean(waterTemp)
TH = LeastSquares3(time(sampleStart));
Tmid = (T0 + TH)/2;

%Find T2 by interpolating the second regression line
tMax = find((LeastSquares2(time) > Tmid), 1, 'first'); %Find index where Tmid intersects the second regression
T2 = LeastSquares3(time(tMax))

startLine = xline(time(sampleStart), 'r--');
midLine = xline(time(tMax), 'b--');
initTemp = plot(time, T0*ones(length(time),1));
hotTemp = plot(time, TH*ones(length(time),1));
midTemp = plot(time, Tmid*ones(length(time),1));
endTemp = plot(time, T2*ones(length(time),1));

subset = [firstReg, secondReg, lastReg, startLine, midLine, initTemp, midTemp, hotTemp, endTemp];

legend(subset, 'First regression', 'Second regression', 'Third regression', "Sample added", "Midpoint temp time", "T0", "Tmid", "TH", "T2", 'Location','best');
%legend("First section","Second section","Third section","Curved section","LS1","LS2","LS3","sampleStart","tMax","T0","TH","Tmid","T2",'Location','best')
ylim([23 30])

%% Best estimate and uncertainty calculations

g = waterTemp;
Dcol1 = time;
Dcol2 = ones([length(time) 1]);
D = [Dcol1 Dcol2];

%Calculate best estimate for the specific heat of sample B
Cs = (Mc*Cc*(T2-T0))/(Ms*(T1-T2))

%Sigma_y based on extrapolation error for each temperature
sigma_T0 = norm([sigma_thermo std(LeastSquares1(firstSection))]);
sigma_T1 = norm([sigma_thermo std(waterTemp)]);
sigma_T2 = norm([sigma_thermo std(LeastSquares3(lastSection))]);

%T0 extrapolation uncertainty
W0 = diag(1./sigma_T0.^2);
extrapT0 = extrapError(W0, A, sampleStart);
sig_Y_T0 = norm([sigma_thermo extrapT0]);

%T1 extrapolation uncertainty
W1 = diag(1./sigma_T1.^2);
extrapT1 = extrapError(W1, D, sampleStart);
sig_Y_T1 = norm([sigma_thermo extrapT1]);

%T2 extrapolation uncertainty
W2 = diag(1./sigma_T2.^2);
extrapT2 = extrapError(W2, C, tMax);
sig_Y_T2 = norm([sigma_thermo extrapT2]);

%Sigma_y based on least squares error for each temperature
LST0 = LSError(LeastSquares1, d, firstSection);
sigma_T0 = norm([sigma_thermo LST0]);

sigma_T1 = norm([sigma_thermo std(waterTemp)]);

LST2 = LSError(LeastSquares3, f, lastSection);
sigma_T2 = norm([sigma_thermo LST2]);

%Choose which sigma_y is larger, least squares error or extrapolated error
if sig_Y_T0 > sigma_T0
    sigma_T0 = sig_Y_T0;
    fprintf("Chose extrapolated error for T0!\n");
else
    fprintf("Chose least squares error for T0!\n");
end

if sig_Y_T1 > sigma_T1
    sigma_T1 = sig_Y_T1;
    fprintf("Chose extrapolated error for T1!\n");
else
    fprintf("Chose least squares error for T1!\n");
end

if sig_Y_T2 > sigma_T2
    sigma_T2 = sig_Y_T2;
    fprintf("Chose extrapolated error for T2!\n");
else
    fprintf("Chose least squares error for T2!\n");
end

%Display all variable uncertainties
sigma_Mc
sigma_Ms
sigma_T0
sigma_T1
sigma_T2

sigma_Cs = calcUncertainty(Cc, Cs, Mc, sigma_Mc, Ms, sigma_Ms, T0, sigma_T0, T1, sigma_T1, T2, sigma_T2)

%sigma_Cs = ((Cs*Cc)/(Ms*(T1-T2)))*sqrt(((T2-T0)*sigma_Mc)^2+(-Mc*(T2-T0)*sigma_Mc/(Ms))^2+(-Mc*sigma_T0)^2+(-Mc*(T2-T0)*sigma_T1/z(T1-T2))^2+(Mc*(T1-T0)*sigma_T2/(T1-T2))^2)

%Display the range of specific heats for sample B
CsHigh = Cs + sigma_Cs
CsLow = Cs - sigma_Cs

%% Percent error calculations

%Calculate percent errors between the best estimate and each given
%material...
percentErrorAc = abs((1.470 - Cs) / 1.470) %if it is actually Acrylic
percentErrorZn = abs((0.402 - Cs) / 0.402) %if it is actually Zinc
percentErrorAl = abs((0.895 - Cs) / 0.895) %if it is actually Aluminum
percentErrorCu = abs((0.261 - Cs) / 0.261) %if it is actually Copper
percentErrorPb = abs((0.129 - Cs) / 0.129) %if it is actually Lead

%Find smallest percent error
percentError = min([percentErrorAc, percentErrorZn, percentErrorAl, percentErrorCu, percentErrorPb])

%Choose the material with the lowest specific heat as the most likely
%sample material
switch percentError
    case percentErrorAc
        fprintf("Smallest percent error is %f, sample is likely Acrylic!\n", percentErrorAc)
    case percentErrorZn
        fprintf("Smallest percent error is %f, sample is likely Zinc!\n", percentErrorZn)
    case percentErrorAl
        fprintf("Smallest percent error is %f, sample is likely Aluminum!\n", percentErrorAl)
    case percentErrorCu
        fprintf("Smallest percent error is %f, sample is likely Copper!\n", percentErrorCu)
    case percentErrorPb
        fprintf("Smallest percent error is %f, sample is likely Lead!\n", percentErrorPb)
    otherwise
        disp("I don't know what this sample is!!!")
end

