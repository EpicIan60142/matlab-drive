clc; clear; close all;

R = 13*0.0254; % in to m
L = 6*0.0254; % in to m
g = 9.81; % m/s^2

measTimes = [5.69, 5.28, 4.94, 2.75]; % sec, trial 2-4, 1
rates = [19, 18, 17, 11] * ((5280*.3048)/3600)/R; % mph to rad/s, trial 2-4, 1

wVec = linspace(min(rates),max(rates), 100);

model = @(w) (2*g*L)./(w*R^2);

figure
hold on
title("Precession Rate vs. Wheel Spin Rate")
plot(rates, (2*pi)./measTimes)
plot(wVec, model(wVec))
xlabel("Wheel spin rate (rad/s)")
ylabel("Precession rate (rad/s)")

legend("Measured Data", "Model Prediction", 'Location', 'best')

figure
hold on
title("Precession Period vs. Wheel Spin Rate")
plot(rates, measTimes)
plot(wVec, (2*pi)./model(wVec))
xlabel("Wheel spin rate (rad/s)")
ylabel("Precession period (sec)")

legend("Measured Data", "Model Prediction", 'Location', 'best')