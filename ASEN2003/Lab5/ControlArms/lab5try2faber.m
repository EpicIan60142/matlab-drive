clear; clc; close all;

syms zeta omega_n Kp Kd real positive

%% Constants

% Aaron found Kp = 5, Kd = 0.5

Mp = 0.2;
settle = 0.05;
ts = 1; % 1 second settling time
step = 5;
%Kp = 1.34; % dervtive gain 
Kg = 33.3; % gear ratio
Km = 0.0401; % poprtional motor constant relating the speed to the motor voltage 
J  = 0.0005+(0.2*0.2794^2)+0.0015; % moment of intertia 
Rm = 19.2; % output resitance of the motor 
%Kd = -0.805; % poprtinal gain 


%% Calculate Gains based on equations
omegaEq = sqrt((Kp*Kg*Km)/(J*Rm));
zetaEq = (Kg^2*Km^2 + Kg*Km*Kd)/(2*sqrt(Kp*Kg*Km*J*Rm));
overshootEq = exp((-1*zeta*pi)/sqrt(1-zeta^2));
settleEq = 1 - (1/sqrt(1-zeta^2))*exp(-zeta*omega_n*ts)*sin(omega_n*sqrt(1-zeta^2)*ts + acos(zeta));

zeta = vpasolve(Mp == exp((-zeta*pi)/sqrt(1-zeta^2)), zeta);
omega_n = solve(1 + settle == 1 - (1/sqrt(1-zeta^2))*exp(-zeta*omega_n*ts)*sin((omega_n*sqrt(1-zeta^2)*ts) + acos(zeta)), omega_n);

Kp = double(vpasolve(omega_n == sqrt((Kp*Kg*Km)/(J*Rm)), Kp));
Kd = double(vpasolve(zeta == (Kg^2*Km^2 + Kg*Km*Kd)/(2*sqrt(Kp*Kg*Km*J*Rm))));

% Kp = 5;
% Kd = 0.5;

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
response = plot(tStep,xStep);
input = plot(tInput, stepInput(tInput, step), '--');
overshoot = yline(step*(1+Mp));
topSettle = yline(step*(1+settle), '--');
bottomSettle = yline(step*(1-settle), '--');
settleLine = xline(ts);

title('Bode plot for time vs angular position')
xlabel('Time (sec)')
ylabel('Angular Position (rad)')
ylim([0, 1 + (1+Mp)*step])

responseLabel = sprintf("Step response: %d%% overshoot, %d%% settling time in %.3f sec\nKp = %.3f, Kd = %.3f", 100*Mp, 100*settle, ts, Kp, Kd);
inputLabel = sprintf("Step input of %.3f", step);

legend([response, input], responseLabel, inputLabel, 'Location', 'best')

