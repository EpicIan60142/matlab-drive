%% ECEN 4638 Closed Loop Post-Processing Script
%
% Ian Faber, Brennen Billig, Luke Hanley
%

%% Processing
close all;

load("VertMotorData.mat")

t = 0:Ts:Te-Ts; % Defome experiment window
t = t + Tt; % shift by Tt

f = 0:1/Te: 1/Ts - 1/Te;
%u = 2*cos(8*pi*t); % Define signal in time

w = 2*pi*f;
% w = w+wr;
% w = [1:N]*wr;

% u = CL_r;
U = fft(u)/length(u); % Perform fft and normalize by number of points

y = theta;
Y = fft(y)/length(y);

H = Y./U;
% H = H./max(H);

s = tf('s');

% %model checks
w0 = 6.28319;
w1 = 3.351;
w2 = w0^2/w1;%11.781;
zeta = 1.0637;
Q = 1/(2*zeta);
K = db2mag(3.0885);
modelsys = tf([K*w0^2], [1 2*zeta*w0 w0^2]);
% modelsys2 = tf([K], [1/(w1*w2) (1/w1)+(1/w2) 1])
% modelsys2 = K/((s/w1 + 1)*(s/w2 + 1));
% modelsys = tf([K], [1/w0^2 1/(Q*w0) 1]);
[modelMag, modelPhase] = bode(modelsys,w);
[modelMag2, modelPhase2] = bode(modelsys2,w);
modelMag = squeeze(db(modelMag));
modelPhase = squeeze(modelPhase);
% modelMag2 = squeeze(db(modelMag2));
% modelPhase2 = squeeze(modelPhase2);

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
% plot(w,modelMag2)
xline(w0,'r--')
xline(1/(2*Ts),'k--')
xlabel("Frequency [rad/s]")
ylabel("Amplitude")
legend("Experiment", "Model", "Nat Freq.", "Nyquist Freq.")
grid on

subplot(2,1,2)
semilogx(w(w<=N*wr), rad2deg(unwrap(angle(H(w<=N*wr)))))
hold on
title("Phase Bode Plot of H")
plot(w,modelPhase);
% plot(w,modelPhase2);
xline(w0,'r--')
xline(1/(2*Ts),'k--')
xlabel("Frequency [rad/s]")
ylabel("Phase [deg]")
legend("Experiment", "Model", "Nat Freq.", "Nyquist Freq.")
grid on


