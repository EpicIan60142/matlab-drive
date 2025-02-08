% ASEN 2003 Lab 2: Bouncing Ball Experiment Error Analysis
% By: Ian Faber, Justin Travis, Niko Von Unger, and Yun Yang

clc; clear; close all;

% Method 2
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