%% HackRF_recieve test ingesting script
% - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Ingest file
fid = fopen("test_hackrf_recieve"); % Open file
raw = fscanf(fid,"%s"); % Convert contents to a vector of chars
fclose(fid); % Close file

numStream = double(raw); % Convert from char to int16

y = numStream(1:2:end) + 1j*numStream(2:2:end); % Convert stream from I's and Q's to a stream of complex numbers

idx = 1:length(y);

Y = fft(real(y));

% binary = dec2bin(numStream); % Convert to binary
% 
% hex = dec2hex(numStream); % Convert to hex
% 
figure
hold on
grid on
plot(real(y(idx)), imag(y(idx)), 'k.') % Plot first 1000 samples
xlim([-256 256])
ylim([-256 256])

figure
plot(Y)
