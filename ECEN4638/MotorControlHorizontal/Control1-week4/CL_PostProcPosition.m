%% ECEN 4638 Closed Loop Post-Processing Script
%
% Ian Faber, Brennen Billig, Luke Hanley
%

%% Processing
close all;

t = 0:Ts:Te-Ts; % Defome experiment window
t = t + Tt; % shift by Tt

f = 0:1/Te: 1/Ts - 1/Te;
%u = 2*cos(8*pi*t); % Define signal in time

w = 2*pi*f;
% w = w+wr;
% w = [1:N]*wr;

u = CL_r;
U = fft(u)/length(u); % Perform fft and normalize by number of points

y = CL_theta;
Y = fft(y)/length(y);

H = Y./U;
% H = H./max(H);

%model checks
modelsys = tf([1.034], [1/(20.1062)^2 1/(1.943*20.1062) 1]);
w0 = 20.1062;
[modelMag, modelPhase] = bode(modelsys,w);
modelMag = squeeze(db(modelMag));
modelPhase = squeeze(modelPhase);

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
subplot(2,1,1)
semilogx(w(w<=N*wr), db(abs(H(w<=N*wr))))
hold on
title("Amplitude Bode Plot of H")
plot(w,modelMag)
xline(w0,'r--')
% xline(1/(2*Ts),'k--')
xlabel("Frequency [rad/s]")
ylabel("Amplitude")
legend("Experiment", "Model")
grid on

subplot(2,1,2)
semilogx(w(w<=N*wr), rad2deg(unwrap(angle(H(w<=N*wr)))))
hold on
title("Phase Bode Plot of H")
plot(w,modelPhase);
xline(w0,'r--')
% xline(1/(2*Ts),'k--')
xlabel("Frequency [rad/s]")
ylabel("Phase [deg]")
legend("Experiment", "Model")
grid on


