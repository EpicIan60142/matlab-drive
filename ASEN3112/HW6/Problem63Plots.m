%% Housekeeping
clc; clear; close all;

%% Constants
P = 10*10^3; % N
EI = 10^6; % Nm^2
L = 0.5; % m

xPlot = linspace(0, L, 1000);

%% Problem 6.3a: Answer
vB = -((4*L^2 - 5*L + 2)*P)/(4*EI*L)
phiB = -((4*L^2 - 6*L + 3)*P)/(2*EI*L^2)

%% Problem 6.3b - Plot v(x) and v'(x)

v = @(x)( (3*(x/L).^2 - 2*(x/L).^3)*vB + ((x/L).^3 - (x/L).^2)*L*phiB );
vPrime = @(x)( (6*x/L - 6*(x/L).^2)*vB + (3*(x/L).^2 - 2*x/L)*L*phiB );

figure
sgtitle("6.3b: Deflection and Slope across Beam")
subplot(2,1,1)
hold on
grid on
title("v(x) vs. x")
plot(xPlot, v(xPlot))
xlabel("x (m)")
ylabel("Deflection (m)")

subplot(2,1,2)
hold on
grid on
title("v'(x) vs. x")
plot(xPlot, vPrime(xPlot))
xlabel("x (m)")
ylabel("Deflection slope (m/m)")

%% Problem 6.3c - Plot M(x) = EI*v''(x)
kappa = @(x)( (6/L - 12*x/L)*vB + (6*x - 2)*phiB );
M = @(x)(EI*kappa(x));

figure
hold on
grid on
title("6.3c: M(x) vs. x")
plot(xPlot, M(xPlot));
xlabel("x (m)")
ylabel("Internal bending moment (N*m)")

