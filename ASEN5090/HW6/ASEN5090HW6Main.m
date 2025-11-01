%% ASEN 5090 HW 6 Main Script
% By: Ian Faber, 11/01/2025

%% Housekeeping
clc; clear; close all;

%% Problem 1: Sampling and frequency resolution
    % Part a. Create a time vector for a sampling frequency of 50 MHz, with
    % a duration of 10 msec.
Fs = 50e6; % Sample frequency of 50 MHz
Ts = 1/Fs; % Sample time
Tdur = 10e-3; % Sample duration of 10 ms

tSample = 0:Ts:Tdur; % Time vector of samples

    % Part b. Nyquist frequency and frequency resolution
fNyquist = Fs/2;
fResolution = Fs/length(tSample)

%% Problem 2: Oscilloscope and spectrum analyzer plots
    % Part a. Make 40 kHz sine wave w/ 1 V amplitude and plot signal +
    % spectrum
fSin = 40e3;
ampSin = 1;
sig1 = ampSin*sin(2*pi*fSin*tSample);

endIdx = 50000;

% figure;
% hold on; grid on;
% title("Oscilloscope Sine Wave Signal")
% plot(tSample(1:endIdx), sig1(1:endIdx), 'b-')
% xlabel("Time [sec]"); ylabel("Signal [V]")

plot_scope_spectrum_analyzer(2, tSample, sig1, tSample(endIdx), [0, 40e4], [-150, 0]);

%% Problem 3: Square wave
    % Part a. Make 40 kHz square wave w/ 1 V amplitude and plot signal +
    % spectrum
fSquare = 40e3;
ampSquare = 1;
sig2 = ampSquare*0.5*(square(2*pi*fSquare*tSample, 50) + 1);

endIdx = 25000;

% figure;
% hold on; grid on;
% title("Oscilloscope Square Wave Signal")
% plot(tSample(1:endIdx), sig2(1:endIdx), 'b-')
% xlabel("Time [sec]"); ylabel("Signal [V]")

plot_scope_spectrum_analyzer(3, tSample, sig2, tSample(endIdx), [0, 40e5]);

%% Problem 4: PRN Codes - Maximal Length
    % Part a. Make a GPS C/A G1 register w/ chipping rate 1.023 MHz
fMaxCode = 1.023e6; % Problem 6e: Try 2 or 0.5 times chipping rate
maxCode = ones(1,1023);
G1 = ones(1,10);
for k = 1:1023
    out1 = G1(end);
    maxCode(k) = out1;
    newBit1 = xor(G1(3), G1(10));
    G1 = [newBit1, G1(1:9)];
end

tMaxCode = zeros(1, length(maxCode));
for k = 0:length(maxCode)
    tMaxCode(k+1) = (1/fMaxCode)*k;
end

tRepeat = length(maxCode)/fMaxCode

timeIdx = 1;
value = maxCode(timeIdx);
sig3 = zeros(1,length(tSample));
for k = 1:length(tSample)
    if tSample(k) >= tMaxCode(timeIdx)
        value = maxCode(timeIdx);
        timeIdx = mod(timeIdx + 1, length(maxCode));
        if timeIdx == 0
            timeIdx = 1;
            tMaxCode = tMaxCode + tRepeat;
        end
    end
    sig3(k) = value;
end

outIdx = 485:495;
time3Out = tSample(outIdx)
sig3Out = sig3(outIdx)

    % Part b. Plot signal as seen by oscilloscope and spectrum analyzer
endIdx = floor(length(tSample)/10);

plot_scope_spectrum_analyzer(4, tSample, sig3, tSample(endIdx), [0, 3*fMaxCode]);

%% Problem 5: PRN Codes - Gold Code
    % Part a. Make GPS C/A code for PRN 9 with chipping rate 1.023 MHz
fGoldCode = 1.023e6;
PRN = 9;
G1 = ones(1,10);
G2 = ones(1,10);
goldCode = generatePRN(G1, G2, PRN, 1023);

tGoldCode = zeros(1, length(goldCode));
for k = 0:length(goldCode)
    tGoldCode(k+1) = (1/fGoldCode)*k;
end

tRepeat = length(maxCode)/fMaxCode

timeIdx = 1;
value = goldCode(timeIdx);
sig4 = zeros(1,length(tSample));
for k = 1:length(tSample)
    if tSample(k) >= tGoldCode(timeIdx)
        value = goldCode(timeIdx);
        timeIdx = mod(timeIdx + 1, length(goldCode));
        if timeIdx == 0
            timeIdx = 1;
            tGoldCode = tGoldCode + tRepeat;
        end
    end
    sig4(k) = value;
end

outIdx = 485:495;
time4Out = tSample(outIdx)
sig4Out = sig4(outIdx)

    % Plot signal and spectrum
endIdx = floor(length(tSample)/10);

plot_scope_spectrum_analyzer(5, tSample, sig4, tSample(endIdx), [0, 3*fGoldCode]);

%% Problem 6: Direct Spread Spectrum Modulation
    % Part a. Make a 5.115 MHz carrier sine wave and BPSK modulate it with
    % the G1 code
fCarrier = 5.115e6; % part d: try 2 or 0.5 times carrier frequency
sigCarrier = sin(2*pi*fCarrier*tSample);
G1PosNeg = convPRNZeroOne2PosNeg(sig3);
sig5 = sigCarrier.*G1PosNeg;

    % Part b. Plot signal and its spectrum
endIdx = floor(length(tSample)/1000);

plot_scope_spectrum_analyzer(6, tSample, sig5, tSample(endIdx), [0, 2*fCarrier]);

%% Problem 7: Demodulation/Carrier Recovery
    % Part a. Add white noise with standard deviation 1 V to 6a's signal






