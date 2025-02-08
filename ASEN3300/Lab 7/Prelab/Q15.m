%% Prelab 7 Q1.5

%% Housekeeping
clc; clear; close all;

%% Constants
Vpp = 3.3;
VDC = 1.65;
minV = VDC - (Vpp/2);
maxV = VDC + (Vpp/2);
bits = 12;
voltages = (Vpp/2)*sin(0:0.1:2*pi) + VDC;


%% Code
[binDec, ~] = voltage2Bin(minV, maxV, bits, voltages)

figure
hold on
grid on
titleText = sprintf("%.0f-bit Bin vs. Array number for %.1f Vpp sine wave with %.2f V DC offset", bits, Vpp, VDC);
title(titleText);
stem(binDec, '.')
xlabel("Array Number")
ylabel("Bin Number")

range = maxV - minV


%% Function
function [binDec, binBin] = voltage2Bin(minV, maxV, bits, voltages)
% voltage2Bin: Determines the bin number, in both decimal and binary, 
%              a given voltage signal would be placed in by an A/D 
%              converter

range = maxV - minV; 
binSize = range/(2^bits);

binDec = floor((voltages-minV)/binSize);
binBin = dec2bin(binDec);

end
