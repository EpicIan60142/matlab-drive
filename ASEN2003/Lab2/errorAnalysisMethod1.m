%% ASEN 2003 Lab 2: Bouncing Ball Experiment Error Analysis
% By: Ian Faber, Justin Travis, Niko Von Unger, and Yun Yang

clc; clear; close all;

%% Method 1
trial1 = readmatrix('.\Method 1 Data\Trial1Data.txt'); % t, y, x
trial2 = readmatrix('.\Method 1 Data\Trial2Data.txt'); % t, y, x
trial3 = readmatrix('.\Method 1 Data\Trial3Data.txt'); % t, x, y
trial4 = readmatrix('.\Method 1 Data\Trial4Data.txt'); % t, x, y
trial5 = readmatrix('.\Method 1 Data\Trial5Data.txt'); % t, x, y

sigHn = 0.0005;
sigH0 = 0.02;

trial1 = trial1(:,2);
trial2 = trial2(:,2);
trial3 = trial3(:,3);
trial4 = trial4(:,3);
trial5 = trial5(:,3);

n = length(trial1);
e1 = (trial1(n)/trial1(1))^(1/(2*n));
sig1 = sqrt( (((1/(2*n*trial1(1)))*(trial1(n)/trial1(1))^((1-2*n)/(2*n)))*sigHn)^2 + (((-1/(2*n*trial1(n)))*(trial1(1)/trial1(n))^((-1-2*n)/(2*n)))*sigH0)^2);

n = length(trial2);
e2 = (trial2(n)/trial2(1))^(1/(2*n));
sig2 = sqrt( (((1/(2*n*trial2(1)))*(trial2(n)/trial2(1))^((1-2*n)/(2*n)))*sigHn)^2 + (((-1/(2*n*trial2(n)))*(trial2(1)/trial2(n))^((-1-2*n)/(2*n)))*sigH0)^2);

n = length(trial3);
e3 = (trial3(n)/trial3(1))^(1/(2*n));
sig3 = sqrt( (((1/(2*n*trial3(1)))*(trial3(n)/trial3(1))^((1-2*n)/(2*n)))*sigHn)^2 + (((-1/(2*n*trial3(n)))*(trial3(1)/trial3(n))^((-1-2*n)/(2*n)))*sigH0)^2);

n = length(trial4);
e4 = (trial4(n)/trial4(1))^(1/(2*n));
sig4 = sqrt( (((1/(2*n*trial4(1)))*(trial4(n)/trial4(1))^((1-2*n)/(2*n)))*sigHn)^2 + (((-1/(2*n*trial4(n)))*(trial4(1)/trial4(n))^((-1-2*n)/(2*n)))*sigH0)^2);

n = length(trial5);
e5 = (trial5(n)/trial5(1))^(1/(2*n));
sig5 = sqrt( (((1/(2*n*trial5(1)))*(trial5(n)/trial5(1))^((1-2*n)/(2*n)))*sigHn)^2 + (((-1/(2*n*trial5(n)))*(trial5(1)/trial5(n))^((-1-2*n)/(2*n)))*sigH0)^2);

w1 = 1/(sig1^2);
w2 = 1/(sig2^2);
w3 = 1/(sig3^2);
w4 = 1/(sig4^2);
w5 = 1/(sig5^2);


eFinal1 = (w1*e1 + w2*e2 + w3*e3 + w4*e4 + w5*e5)/(w1 + w2 + w3 + w4 + w5)
sigFinal1 = 1/sqrt(w1 + w2 + w3 + w4 + w5)


%% Method 2
method2 = readmatrix('.\Method 2 Data\Delta.xlsx');


%% Method 3
method3 = readmatrix('.\Method 3 Data\Method 3 Data.xlsx');

