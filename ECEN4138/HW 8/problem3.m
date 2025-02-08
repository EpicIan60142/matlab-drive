%% ECEN 4138 HW 8 Problem 5.26
% - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
zeta = 0.607;
theta = (pi/2 + asin(zeta):0.001:(3*pi/2)-asin(zeta));

wn = 4.5; % rad/s

t = -100:0.001:-wn*sqrt(1-zeta^2)*tan(zeta);

% Lead compensator
z1 = 6;
p1 = 5*z1;
K = 700;

% Lag compensator
z2 = 0.24;
p2 = 0.01;

alpha = (-1-10-p1-p2+z1+z2)/3;

%% Define L(s)

s = tf('s');

G = 10/(s*(s+1)*(s+10));
C = ((s+z1)^2*(s+z2))/((s+p1)^2*(s+p2));

% L1 = (10*(s+z1)^2*(s+z2))/(s*(s+1)*(s+10)*(s+p1)^2*(s+p2)); % L(s) for 5.20
L1 = G*C;

L = L1; % Choose root locus

%% Plot root locus
figure
hold on

rlocus(L)

% Rise time requirement
y = wn*sin(theta);
x = wn*cos(theta);
plot(x,y, 'k--');

% Overshoot requirement
z = tan(zeta)^-1*t;
a = plot(t,z, 'k--');
plot(t,-z,'k--')

b = xline(alpha, 'b--');

xlim([-1.5*p1, 10])
ylim([-30, 30])

legend([a,b], ["Performance Region", "Center of Asymptotes"], 'Location', 'best');

%% Simulate responses
C = K*C;
L = G*C;

T = minreal(feedback(L,1))

time = (0:0.001:100)';

figure
step(T); % Simulate step response

figure
resp = lsim(T, time, time); % Simulate ramp response
lsim(T, time, time)

disp(time(end) - resp(end));


