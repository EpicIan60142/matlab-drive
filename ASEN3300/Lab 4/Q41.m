%% ASEN 3300 Lab 4 Section 4.1 script

%% Housekeeping
clc; clear; close all;

%% Code
signal1 = load('lab04_analysis_signal1.mat');
signal2 = load('lab04_analysis_signal2.mat');
signal3 = load("lab04_analysis_signal3.mat");

%% 4.1a
signal = signal1;
nPoints = 100;%length(signal.x);

figure
hold on
title("Signal 1, first 100 points")
plot(signal.t(1:nPoints), signal.x(1:nPoints));
xlabel("Time [sec]")
ylabel("Voltage [V]")

%% 4.1b
% freq = linspace(0, signal1.Fs/2, length(signal1.x(1:100)));
% psd1 = fft(signal1.x(1:100));
x = signal.x(1:nPoints);
fs = signal.Fs;

N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
% psdx = psdx/sqrt(2);
% psdx = psdx.^2;
freq = 0:fs/length(x):fs/2;

figure
hold on
title("Signal 1 Power Spectral Density, first 100 points")
plot(freq,pow2db(psdx))
grid on
xlabel("Frequency [Hz]")
ylabel("Power/Frequency [dB($\frac{Vrms^2}{Hz}$)]",'Interpreter','latex')

%% 4.1c
x = signal.x;
fs = signal.Fs;

N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
% psdx = psdx/sqrt(2);
% psdx = psdx.^2;
freq = 0:fs/length(x):fs/2;

figure
hold on
title("Signal 1 Power Spectral Density, full signal")
plot(freq,pow2db(psdx))
grid on
xlabel("Frequency [Hz]")
ylabel("Power/Frequency [dB($\frac{Vrms^2}{Hz}$)]",'Interpreter','latex')

%% 4.1d
x = signal2.x;
fs = signal2.Fs;

figure
hold on
title("Signal 2, full signal")
plot(signal2.t, signal2.x);
xlabel("Time [sec]")
ylabel("Voltage [V]")

N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
% psdx = psdx/sqrt(2);
% psdx = psdx.^2;
freq = 0:fs/length(x):fs/2;

figure
hold on
title("Signal 2 Power Spectral Density, full signal")
plot(freq,pow2db(psdx))
grid on
xlabel("Frequency [Hz]")
ylabel("Power/Frequency [dB($\frac{Vrms^2}{Hz}$)]",'Interpreter','latex')

%% 4.1e
x = signal3.x;
fs = signal3.Fs;

figure
hold on
title("Signal 3, full signal")
plot(signal3.t, signal3.x);
xlabel("Time [sec]")
ylabel("Voltage [V]")

N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
% psdx = psdx/sqrt(2);
% psdx = psdx.^2;
freq = 0:fs/length(x):fs/2;

figure
hold on
title("Signal 3 Power Spectral Density, full signal")
plot(freq,pow2db(psdx))
grid on
xlabel("Frequency [Hz]")
ylabel("Power/Frequency [dB($\frac{Vrms^2}{Hz}$)]",'Interpreter','latex')




