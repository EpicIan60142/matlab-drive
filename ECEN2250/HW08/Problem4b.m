clear; clc; close all;

I0 = -1;
R = 3000;
w = 10e6;
L = 2/3000;

tau = L/R;

time = 0:tau/1000:10*tau;

current = @(t) (I0 + (2*w*tau)/(1+(w*tau)^2))*exp(-t/tau) + (2/(1+(w*tau)^2))*sin(w*t) - ((2*w*tau)/(1+(w*tau)^2))*cos(w*t);
natural = @(t) (I0 + (2*w*tau)/(1+(w*tau)^2))*exp(-t/tau);
forced = @(t) (2/(1+(w*tau)^2))*sin(w*t) - ((2*w*tau)/(1+(w*tau)^2))*cos(w*t);

figure
hold on
title("First-Order Circuit Response for Sinusoidal Input")
plot(time, current(time), 'b-');
plot(time, natural(time), 'r--');
plot(time, forced(time), 'k--');
xlabel("time (sec)");
ylabel("Current (A)");
legend("Full current response", "Natural current response", "Forced current response", 'Location', 'best');


