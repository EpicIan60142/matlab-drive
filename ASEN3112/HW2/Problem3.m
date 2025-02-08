clc; clear; close all;

delta = (0:1:1000)/1000; % mm to m
theta = acosd(4/5); % Degrees
L0 = 5; % m
Ab = 0.001; % m^2
E = 100*(10^6); % Mpa to Pa

approxP = @(delta)(25600*delta);

L = @(delta)(sqrt(L0^2 + delta.^2 - 2*L0*delta*cosd(180-theta)));
phi = @(L)(asind((L0*sind(180-theta)./L))); % Degrees
trueP = @(L, phi)(2*E*Ab*((L-L0)/L0).*cosd(phi));

pError = @(trueP, approxP)(100*((approxP - trueP)./trueP));

a = approxP(delta);
b = trueP(L(delta), phi(L(delta)));

figure()

subplot(1,2,1)
hold on
grid on
title("Predictions of Force P vs. Displacement")
xlabel("Displacement(m)")
ylabel("Force (N)")
plot(delta, approxP(delta))
plot(delta, trueP(L(delta), phi(L(delta))))
legend("Small Displacement Prediction", "General Prediction")
hold off

subplot(1,2,2)
hold on
grid on
title("Force Prediction Error vs. Displacement")
xlabel("Displacement (m)")
ylabel("Error (% of general prediction)")
plot(delta, pError( trueP(L(delta), phi(L(delta))) , approxP(delta) ) )
hold off


