%% ECEN 4638 Post-Processing Script - Velocity Data
%
% Ian Faber, Brennen Billig, Luke Hanley
%

%% Housekeeping
close all;

%% Load Data (for offline analysis)
load("motorSystem_Data.mat")

%% Processing
t = 0:Ts:Te-Ts; % Define experiment window
t = t + Tt; % shift by Tt

f = 0:1/Te: 1/Ts - 1/Te;
w = 2*pi*f;

u = input
U = fft(u)/length(u); % Perform fft and normalize by number of points

y = outputVel
Y = fft(y)/length(y);

H = Y./U;

% figure
% hold on;
% title("signals")
% plot(t,u,t,y)
% xlabel("Time [sec]")
% ylabel("Amplitude")
% legend("u","y")
% 
% figure
% hold on;
% title("FFT of u")
% stem(f,abs(U))
% xline(1/(2*Ts),'k--')
% xlabel("Frequency [Hz]")
% ylabel("Amplitude")

figure
hold on;
title("FFT of y")
stem(f,abs(Y))
xline(1/(2*Ts),'k--')
xlabel("Frequency [Hz]")
ylabel("Amplitude")

%model checks
modelsys = tf([4.9], [0.0886 1]);
w0 = 1/0.0886;
[modelMag, modelPhase] = bode(modelsys,w);
modelMag = squeeze(db(modelMag));
modelPhase = squeeze(modelPhase);

figure
subplot(2,1,1)
semilogx(w(w<=N*wr), db(abs(H(w<=N*wr))))
hold on
title("Amplitude Bode Plot of H")
plot(w,modelMag)
yline(12.7601-3,'b--')
yline(13.7801-3,'r--')
xline(11.629, 'b--')
xline(w0,'r--')
xlabel("Frequency [rad/s]")
ylabel("Amplitude")
legend("Experiment", "Model")
grid on

subplot(2,1,2)
semilogx(w(w<=N*wr), rad2deg(angle(H(w<=N*wr))))
hold on
title("Phase Bode Plot of H")
plot(w,modelPhase);
xlabel("Frequency [rad/s]")
ylabel("Phase [deg]")
legend("Experiment", "Model")
grid on


