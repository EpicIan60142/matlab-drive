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

trial1 = method2(1:3,1);
trial2 = method2(1:3,2);
trial3 = method2(1:3,3);
trial4 = method2(1:3,4);
trial5 = method2(1:3,5);

n=length(trial1);

e1 = (trial1(n)/trial1(n-1));
e2 = (trial2(n)/trial2(n-1));
e3 = (trial3(n)/trial3(n-1));
e4 = (trial4(n)/trial4(n-1));
e5 = (trial5(n)/trial5(n-1));


sig1 = .01;

sig_bounce1 = sqrt(((sig1/(trial1(n-1))^2 + ((-trial1(n)*sig1)/(trial1(n-1))^2))^2));
sig_bounce2 = sqrt(((sig1/(trial2(n-1))^2 + ((-trial2(n)*sig1)/(trial2(n-1))^2))^2));
sig_bounce3 = sqrt(((sig1/(trial3(n-1))^2 + ((-trial3(n)*sig1)/(trial3(n-1))^2))^2));
sig_bounce4 = sqrt(((sig1/(trial4(n-1))^2 + ((-trial4(n)*sig1)/(trial4(n-1))^2))^2));
sig_bounce5 = sqrt(((sig1/(trial5(n-1))^2 + ((-trial5(n)*sig1)/(trial5(n-1))^2))^2));

w1 = 1/(sig_bounce1^2);
w2 = 1/(sig_bounce2^2);
w3 = 1/(sig_bounce3^2);
w4 = 1/(sig_bounce4^2);
w5 = 1/(sig_bounce5^2);

eFinal2 = (w1*e1 + w2*e2 + w3*e3 + w4*e4 + w5*e5)/(w1 + w2 + w3 + w4 + w5)
sigFinal2 = 1/sqrt(w1 + w2 + w3 + w4 + w5)


%% Method 3
method3 = readmatrix('.\Method 3 Data\Method 3 Data.xlsx');

Trial1 = method3(1,2);
Trial2 = method3(2,2);
Trial3 = method3(3,2);
Trial4 = method3(4,2);
Trial5 = method3(5,2);

ts1 = Trial1;
ts2 = Trial2;
ts3 = Trial3;
ts4 = Trial4;
ts5 = Trial5;

g = 9.81;
n = length(method3(1:5,2));
h_0 = 1; %meter
sig_ts = 0.01;
sig_h_0 = 0.02;


Sig_stop_Trial1 = sqrt((2*sqrt(2*h_0/g)*(sig_ts)/(ts1 + sqrt(2*h_0/g)^2)^2)+((2*ts1*(sig_h_0)^2)/(g*sqrt(2*h_0/g)*(ts1+sqrt(2*h_0/g)^2))));
Sig_stop_Trial2 = sqrt((2*sqrt(2*h_0/g)*(sig_ts)/(ts2 + sqrt(2*h_0/g)^2)^2)+((2*ts2*(sig_h_0)^2)/(g*sqrt(2*h_0/g)*(ts2+sqrt(2*h_0/g)^2))));
Sig_stop_Trial3 = sqrt((2*sqrt(2*h_0/g)*(sig_ts)/(ts3 + sqrt(2*h_0/g)^2)^2)+((2*ts3*(sig_h_0)^2)/(g*sqrt(2*h_0/g)*(ts3+sqrt(2*h_0/g)^2))));
Sig_stop_Trial4 = sqrt((2*sqrt(2*h_0/g)*(sig_ts)/(ts4 + sqrt(2*h_0/g)^2)^2)+((2*ts4*(sig_h_0)^2)/(g*sqrt(2*h_0/g)*(ts4+sqrt(2*h_0/g)^2))));
Sig_stop_Trial5 = sqrt((2*sqrt(2*h_0/g)*(sig_ts)/(ts5 + sqrt(2*h_0/g)^2)^2)+((2*ts5*(sig_h_0)^2)/(g*sqrt(2*h_0/g)*(ts5+sqrt(2*h_0/g)^2))));

e1 = ((ts1 - sqrt(2*h_0/g))/(ts1 + sqrt(2*h_0/g)));
e2 = ((ts2 - sqrt(2*h_0/g))/(ts2 + sqrt(2*h_0/g)));
e3 = ((ts3 - sqrt(2*h_0/g))/(ts3 + sqrt(2*h_0/g)));
e4 = ((ts4 - sqrt(2*h_0/g))/(ts4 + sqrt(2*h_0/g)));
e5 = ((ts5 - sqrt(2*h_0/g))/(ts5 + sqrt(2*h_0/g)));

w1 = 1/((Sig_stop_Trial1)^2);
w2 = 1/((Sig_stop_Trial2)^2);
w3 = 1/((Sig_stop_Trial3)^2);
w4 = 1/((Sig_stop_Trial4)^2);
w5 = 1/((Sig_stop_Trial5)^2);

eFinal3 = (w1*e1 + w2*e2 + w3*e3 + w4*e4 + w5*e5)/(w1 + w2 + w3 + w4 + w5)
sigFinal3 = 1/sqrt(w1 + w2 + w3 + w4 + w5)

%% Combined, final e

w1 = 1/(sigFinal1)^2;
w2 = 1/(sigFinal2)^2;
w3 = 1/(sigFinal3)^2;

eFinal = (w1*eFinal1 + w2*eFinal2 + w3*eFinal3)/(w1 + w2 + w3)
sigFinal = 1/sqrt(w1 + w2 + w3)

%% Plotting

coefficients = [eFinal1; eFinal2; eFinal3];
sigmas = [sigFinal1; sigFinal2; sigFinal3];

figure
errorbar(coefficients, sigmas, 'ro')

xlim([0.5 3.5])
ylim([0.775 0.8125])
xlabel("Method")
ylabel("Coefficient of restitution")
title("Coefficients of Restitution Using Different Dropping Methods")

