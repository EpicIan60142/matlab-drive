%% ECEN 4638 Post-Processing Script - Position Data
%
% Ian Faber, Brennen Billig, Luke Hanley
%

%% Housekeeping
close all;

%% Load Data (for offline analysis)
load("flexJointSysID(no ian).mat")

%% Processing
t = 0:Ts:Te-Ts; % Define experiment window
t = t + Tt; % shift by Tt

f = 0:1/Te: 1/Ts - 1/Te;
w = 2*pi*f;

% ------------ H1: U to Theta 2 -------------- %

% volt = u;
u1 = volt;
U1 = fft(u1)/length(u1); % Perform fft and normalize by number of points

y1 = theta2;
Y1 = fft(y1)/length(y1);

H1 = Y1./U1; % u to theta2 -> theta2/u

% ------------ H2: Theta 1 to Theta 2 -------- %

u2 = theta1;
U2 = fft(u2)/length(u2);

y2 = theta2;
Y2 = fft(y2)/length(y2);

H2 = Y2./U2; % theta1 to theta2 -> theta2/theta1



% --------------- H3: U to Theta 1  --------- %

% To get u to theta1, need theta2/u * theta1/theta2
u3 = volt;
U3 = fft(u3)/length(u3);

y3 = theta1;
Y3 = fft(y3)/length(y3);


H3 = Y3./U3;

% working with H1, u to theta 1
% H3 = H1;            

%zero at origin, subtract off H1;
% s = tf('s');
% ogz = s;
% H3 = H3 - ogz;


H4 = H1./H2; % Check if SysID is right

% figure
% hold on;
% title("signals")
% plot(t,u1,t,y1)
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

% figure
% hold on;
% title("FFT of y")
% stem(f,abs(Y))
% xline(1/(2*Ts),'k--')
% xlabel("Frequency [Hz]")
% ylabel("Amplitude")

%model checks
% K = db2mag(6.023);
% w0 = 3.351;
% Q = db2mag(-37.4798 - K);
% % modelsys = tf([K], [1/(w0^2) 1/(Q*w0) 1]);
% % modelsys = tf(1,[1 0])
% modelsys = tf([K], [1/w0 1 0]);
% % w0 = 28.9027;
% [modelMag, modelPhase] = bode(modelsys,w);
% modelMag = squeeze(db(modelMag));
% modelPhase = squeeze(modelPhase);

% wt = 16.3363;

% Black Boxing H1
s = tf('s');
% model1 = (474.4 * s)/((s + 37.07)*(s - 37.07) * (s + 52.15));
% model1 = (0.0065 * s);                % just zero at OG
% model1 = (9.08 * s)/(s^2 + 22.06*s + 37.07^2);  
model1 = (154 * s)/((s^2 + 20.83*s + 35^2)*(s + 18));


[model1Mag, model1Phase] = bode(model1,w);
model1Mag = squeeze(db(model1Mag));
model1Phase = squeeze(model1Phase);

% Black Boxing H2
model2 = (1 * s^2)/(s^2 + 1.71 * s + 25.76^2);

[model2Mag, model2Phase] = bode(model2,w);
model2Mag = squeeze(db(model2Mag));
model2Phase = squeeze(model2Phase);

% Black Boxing H3

model3 = model1 / model2;
[model3Mag, model3Phase] = bode(model3,w);
model3Mag = squeeze(db(model3Mag));
model3Phase = squeeze(model3Phase);

% ---------------------- H1 (u to theta 2) --------------- %
figure
subplot(2,1,1)
semilogx(w(w<=N*wr), db(abs(H1(w<=N*wr)))) % Experiment
% semilogx(w(w<=N*wr), db(abs(H4(w<=N*wr)))) % Experiment
hold on
title("Amplitude Bode Plot of H1 (u to Theta 2)")
plot(w,model1Mag ) % Model
% yline(13.17-db(wt)-3,'b--') % Experiment w0
% yline(-3.91276-db(w0)-3,'r--') % Model w0
% xline(wt, 'b--') % Experiment w0
% xline(w0,'r--') % Model w0
xlabel("Frequency [rad/s]")
ylabel("Amplitude [dB]")
legend("Raw H1", "model")
grid on

subplot(2,1,2)
semilogx(w(w<=N*wr), rad2deg(unwrap(angle(H1(w<=N*wr))))) % Experiment
% semilogx(w(w<=N*wr), rad2deg(unwrap(angle(H4(w<=N*wr))))) % Experiment
hold on
title("Phase Bode Plot of H1 (u to Theta 2)")
plot(w,model1Phase); % Model
% xline(w0,'r--') % Model w0
xlabel("Frequency [rad/s]")
ylabel("Phase [deg]")
legend("Raw H1", "model")
grid on

% ---------------------- H2 (theta 1 to theta 2) --------------- %
figure
subplot(2,1,1)
semilogx(w(w<=N*wr), db(abs(H2(w<=N*wr)))) % Experiment
% semilogx(w(w<=N*wr), db(abs(H4(w<=N*wr)))) % Experiment
hold on
title("Amplitude Bode Plot of H2 (Theta 1 to Theta 2)")
plot(w,model2Mag) % Model
% yline(13.17-db(wt)-3,'b--') % Experiment w0
% yline(-3.91276-db(w0)-3,'r--') % Model w0
% xline(wt, 'b--') % Experiment w0
% xline(w0,'r--') % Model w0
xlabel("Frequency [rad/s]")
ylabel("Amplitude [dB]")
legend("Raw H2", "model")
grid on

subplot(2,1,2)
semilogx(w(w<=N*wr), rad2deg(unwrap(angle(H2(w<=N*wr)))) + 360) % Experiment
% semilogx(w(w<=N*wr), rad2deg(unwrap(angle(H4(w<=N*wr))))) % Experiment
hold on
title("Phase Bode Plot of H2 (Theta 1 to Theta 2)")
plot(w,model2Phase); % Model
% xline(w0,'r--') % Model w0
xlabel("Frequency [rad/s]")
ylabel("Phase [deg]")
legend("Raw H2", "model")
grid on


% ------------------------ H3 (u to theta 1) --------------%

figure
subplot(2,1,1)
semilogx(w(w<=N*wr), db(abs(H3(w<=N*wr)))) % Experiment
% semilogx(w(w<=N*wr), db(abs(H4(w<=N*wr)))) % Experiment
hold on
title("Amplitude Bode Plot of H3 (u to Theta 1)")
plot(w,model3Mag) % Model
% yline(13.17-db(wt)-3,'b--') % Experiment w0
% yline(-3.91276-db(w0)-3,'r--') % Model w0
% xline(wt, 'b--') % Experiment w0
% xline(w0,'r--') % Model w0
xlabel("Frequency [rad/s]")
ylabel("Amplitude [dB]")
legend("Raw H3", "model")
grid on

subplot(2,1,2)
semilogx(w(w<=N*wr), rad2deg(unwrap(angle(H3(w<=N*wr)))) - 360) % Experiment
% semilogx(w(w<=N*wr), rad2deg(unwrap(angle(H4(w<=N*wr))))) % Experiment
hold on
title("Phase Bode Plot of H3 (u to Theta 1)")
plot(w,model3Phase); % Model
% xline(w0,'r--') % Model w0
xlabel("Frequency [rad/s]")
ylabel("Phase [deg]")
legend("Raw H3", "model")
grid on


% ------------- State Space ----------- %

T = [model3; model1; s.*model3; s.*model1;];
oldSys = ss(T,'minimal')

A = oldSys.A;
B = oldSys.B;
C = oldSys.C;
D = oldSys.D;

Ahat = C*A*C^-1;
Bhat = C*B;
Chat = eye(4);
Dhat = zeros(4,1);

newSys = ss(Ahat, Bhat, Chat, Dhat) % Open Loop TF


