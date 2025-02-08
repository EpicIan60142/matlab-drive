%% ECEN 4638 Post-Processing Script - Position Data
%
% Ian Faber, Brennen Billig, Luke Hanley
%

%% Housekeeping
close all;

%% Load Data (for offline analysis)
load("freqDataNoSpin.mat")
volt = u; % Fix voltage naming

%% Processing
t = 0:Ts:Te-Ts; % Define experiment window
t = t + Tt; % shift by Tt

f = 0:1/Te: 1/Ts - 1/Te;
w = 2*pi*f;

% volt = u;
u = volt;
U = fft(u)/length(u); % Perform fft and normalize by number of points

y = theta;
Y = fft(y)/length(y);

H = Y./U;

figure
hold on;
title("signals")
plot(t,u,t,y)
xlabel("Time [sec]")
ylabel("Amplitude")
legend("u","y")
% 
% figure
% hold on;
% title("FFT of u")
% stem(f,abs(U))
% xline(1/(2*Ts),'k--')
% xlabel("Frequency [Hz]")
% ylabel("Amplitude")

% figure
% hold on;
% title("FFT of y")
% stem(f,abs(Y))
% xline(1/(2*Ts),'k--')
% xlabel("Frequency [Hz]")
% ylabel("Amplitude")

%model checks
K = db2mag(6.023)
w0 = 3.351;
Q = db2mag(-37.4798 - K);
% modelsys = tf([K], [1/(w0^2) 1/(Q*w0) 1]);
% modelsys = tf(1,[1 0])
modelsys = tf([K], [1/w0 1 0]);
% w0 = 28.9027;
[modelMag, modelPhase] = bode(modelsys,w);
modelMag = squeeze(db(modelMag));
modelPhase = squeeze(modelPhase);

% wt = 16.3363;

figure
subplot(2,1,1)
semilogx(w(w<=N*wr), db(abs(H(w<=N*wr)))) % Experiment
hold on
title("Amplitude Bode Plot of H")
plot(w,modelMag) % Model
% yline(13.17-db(wt)-3,'b--') % Experiment w0
% yline(-3.91276-db(w0)-3,'r--') % Model w0
% xline(wt, 'b--') % Experiment w0
% xline(w0,'r--') % Model w0
xlabel("Frequency [rad/s]")
ylabel("Amplitude [dB]")
legend("Experiment", "Model")
grid on

subplot(2,1,2)
semilogx(w(w<=N*wr), rad2deg(unwrap(angle(H(w<=N*wr))))) % Experiment
hold on
title("Phase Bode Plot of H")
plot(w,modelPhase); % Model
% xline(w0,'r--') % Model w0
xlabel("Frequency [rad/s]")
ylabel("Phase [deg]")
legend("Experiment", "Model")
grid on


