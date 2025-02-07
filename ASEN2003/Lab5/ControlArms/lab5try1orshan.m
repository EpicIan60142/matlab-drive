clear; clc; close all;

syms zeta omega_n Kp Kd real positive

%% Constants

% Aaron found Kp = 5, Kd = 0.5

step = 0.5;
Kg = 33.3; % gear ratio
Km = 0.0401; % poprtional motor constant relating the speed to the motor voltage 
J  = 0.0005+(0.2*0.2794^2)+0.0015; % moment of intertia 
Rm = 19.2; % output resitance of the motor 

Kp = 5;
Kd = 0.5;

%% Closed Loop System
d2 = 1;                                        % term connected to s^2
d1 = (Kg^2 * Km^2)/(J*Rm) + (Kd*Kg*Km)/(J*Rm); % term connected to s
d0 = (Kp*Kg*Km)/(J*Rm);                        % term connected to 1

num = (Kp*Kg*Km)/(J*Rm);
den = [d2 d1 d0];
sysTF = tf(num,den);

sinInput = @(t) sin(t);
rampInput = @(t) t;
stepInput = @(t,c) c*ones(length(t),1);
tInput = 0:1/1000:5;

%% Step Response
[xSin,tSin] = lsim(sysTF, sinInput(tInput), tInput);
[xRamp, tRamp] = lsim(sysTF, rampInput(tInput), tInput);
[xStep, tStep] = lsim(sysTF, stepInput(tInput, step), tInput);


figure(1)
hold on
plot(tStep,xStep)
plot(tInput, stepInput(tInput, step), '--')
yline(step*(1+Mp));
yline(step*(1+settle), '--')
yline(step*(1-settle), '--')
xline(1)
title('Bode plot for time vs angular position')
xlabel('Time (sec)')
ylabel('Angular Position (rad)')
legend("response", "input")

