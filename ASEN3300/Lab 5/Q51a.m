%% ASEN 3300 Lab 5 Q5.1a

%% Housekeeping
clc; clear; close all

%% Code

Vin = [1,2,3,4,5,5.6,6,7,8,9,10];
Vdiode = [0.989,1.992,2.995,3.989,4.953,5.332,5.431,5.505,5.536,5.553,5.571];

figure
hold on
title("Zener Diode Output vs. Input Voltage")
plot(Vin, Vdiode, '.', 'MarkerSize',15);
yline(5.6,'k--')
xlabel("Input Voltage [V]")
ylabel("Output Voltage [V]")

legend("Voltage","Theoretical Diode Voltage Drop",'location','best')

