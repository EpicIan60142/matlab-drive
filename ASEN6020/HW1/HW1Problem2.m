%% ASEN 6020 HW 1 Problem 2 Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Cost equations
J_BE = @(r, l) sqrt((2.*l)./(1+l)) - 1 + sqrt((2.*r)./(l.*(l+r))) - sqrt(2./(l.*(l+1))) + (1./sqrt(r)) - sqrt((2.*l)./(r.*(l+r))); 

J_H = @(r) sqrt((2.*r)./(1+r)) - 1 + 1./sqrt(r) - sqrt(2./(r.*(1+r)));

%% Define r and l arrays
r = 1:0.1:10;
l = 1:0.1:10;

[R, L] = meshgrid(r,l);

mask = L<=R;

cost = J_BE(R(mask),L(mask)) - J_H(R(mask));

figure
hold on; grid on;
title("Hohmann Transfer Cost vs. r")
plot(r,J_H(r));
xlabel("r"); ylabel("J_H")

figure
hold on; grid on;
title("Transfer Cost Comparison vs. r, l")
plot3(R(mask),L(mask),J_BE(R(mask),L(mask)), '.')
plot3(R(mask),L(mask),J_H(R(mask)),'.')
xlabel("R"); ylabel("L"); zlabel("J_{BE}")
legend("J_{BE}", "J_H"); view([30 35])

figure
hold on; grid on;
title("Difference in Cost Between J_{BE} and J_H (J_{BE} - J_H)")
plot3(R(mask),L(mask),cost,'.')
xlabel("R"); ylabel("L"); zlabel("cost")
view([30 35])

