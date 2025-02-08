%% ECEN 4638 Post-Processing Script
%
% Ian Faber, Brennen Billig, Luke Hanley
%

%% Processing
close all;

load("frequencyData.mat")

t = 0:Ts:T-Ts; % Defome experiment window
t = t + Tt; % shift by Tt

f = 0:1/T: 1/Ts - 1/T;
%u = 2*cos(8*pi*t); % Define signal in time

w = 2*pi*f;
% w = w+wr;
% w = [1:N]*wr;

u = out.u.Data
U = fft(u)/length(u); % Perform fft and normalize by number of points

y = out.y.Data
Y = fft(y)/length(y);

H = Y./U;
% H = H./max(H);

figure
hold on;
title("signals")
plot(t,u,t,y)
xlabel("Time [sec]")
ylabel("Amplitude")
legend("u","y")

figure
hold on;
title("FFT of u")
stem(f,abs(U))
xline(1/(2*Ts),'k--')
xlabel("Frequency [Hz]")
ylabel("Amplitude")

figure
hold on;
title("FFT of y")
stem(f,abs(Y))
xline(1/(2*Ts),'k--')
xlabel("Frequency [Hz]")
ylabel("Amplitude")

figure
title("Amplitude Bode Plot of H")
semilogx(w(w<=N*wr), db(abs(H(w<=N*wr))))
% xline(1/(2*Ts),'k--')
xlabel("Frequency [rad/s]")
ylabel("Amplitude")
hold on
grid on

figure
title("Phase Bode Plot of H")
semilogx(w(w<=N*wr), rad2deg(angle(H(w<=N*wr))))
% xline(1/(2*Ts),'k--')
xlabel("Frequency [rad/s]")
ylabel("Phase [deg]")
grid on
