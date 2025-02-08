%% ASEN 4013 HW 10 Problem 1
% - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
a = 0.012; % in/s
n = 0.45; % n.d.

% Pc = linspace(0.101, 161.1, 1000);

% r = a*Pc.^n;
r = 0.139658; % in/s
% r = 0.15

% t = linspace(0,3,length(Pc));
t = 0:0.001:3;
t1 = 0.299487./r;

sigma1 = 0.435*(1-cos(asin((0.016+r.*t)/0.435)));
sigma2 = sqrt(0.435^2-r.^2.*t.^2);

sigma3 = acos((0.016+r.*t)/0.435);
sigma4 = asin(sigma2/0.435);

sigma5 = asin((0.016+r.*t)/0.435);

%% Calculate areas
Aslot = (t<t1).*(1.72733 + 6.758*r.*t - 4*r.^2.*t.^2 + sigma1.*(2*r.*t - 3.83) ) + ...
        (t>=t1).*(3.3321 - 1.740*r.*t + sigma1.*(4*r.*t - 7.66) + sigma2.*(3.83 - 2*r.*t));

figure
hold on
title("A_{slot} vs. time")
plot(t, Aslot)
xlabel("Time [sec]")
ylabel("Area [in^2]")

Aaft = (t<t1).*(0.435^2*(pi - sigma5) - 2*r.^2.*t.^2 - (0.032 + 0.435 - sigma1).*r.*t - 0.016*0.435 + 0.016*sigma1) + ...
       (t>=t1).*(0.435^2*(2*sigma3 + sigma4) - 0.87*r.*t - 0.01392 + sigma1.*(2*r.*t + 0.032) - sigma2.*r.*t);

figure
hold on
title("A_{aft} vs. time")
plot(t, Aaft)
xlabel("Time [sec]")
ylabel("Area [in^2]")

feasible = Aaft >= 0 & Aslot >= 0; % If either area term is 0, we've burned through all the propellant.

Ab = Aslot(feasible) + Aaft(feasible);
t = t(feasible);

figure
hold on
title("A_b vs. time")
plot(t, Ab)
xlabel("Time [sec]")
ylabel("Area [in^2]")

%% Calculate thrust and impulse
F = 2.169*Ab;

peakThrust = max(F)
avgThrust = mean(F)

totalImpulse = trapz(t, F) % lbf-sec

figure
hold on
title("Thrust vs. time")
plot(t,F)
xlabel("Time [sec]")
ylabel("Thrust [lb_f]")



