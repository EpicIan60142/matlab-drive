%% Prelab 7 Q1.4

%% Housekeeping
clc; clear; close all;

%% Constants
minV = 0;
maxV = 3.3;
voltages = 0:0.25:3.25;

%% 4 bits
bits = 4;
[binDec4, binBin4] = voltage2Bin(minV, maxV, bits, voltages);

figure
hold on
grid on
titleText = sprintf("Bin number vs. Voltage for %.0f-bit ADC", bits);
title(titleText);
stem(voltages, binDec4, '-')
xlabel("Voltages [V]")
ylabel("Bin Number")

%% 8 bits
bits = 8;
[binDec8, binBin8] = voltage2Bin(minV, maxV, bits, voltages);

figure
hold on
grid on
titleText = sprintf("Bin number vs. Voltage for %.0f-bit ADC", bits);
title(titleText);
stem(voltages, binDec8, '.')
xlabel("Voltages [V]")
ylabel("Bin Number")

%% 12 bits
bits = 12;
[binDec12, binBin12] = voltage2Bin(minV, maxV, bits, voltages);

figure
hold on
grid on
titleText = sprintf("Bin number vs. Voltage for %.0f-bit ADC", bits);
title(titleText);
stem(voltages, binDec12, '.')
xlabel("Voltages [V]")
ylabel("Bin Number")

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
