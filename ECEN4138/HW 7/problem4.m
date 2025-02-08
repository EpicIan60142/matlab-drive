%% ECEN 4138 HW 7 Problem 5.13
% - Ian Faber

%% Housekeeping
clc; clear; close all;

addpath('..\')

%% Define L(s)

s = tf('s');

L1 = ((s+1)*(s^2+81))/(s^2*(s^2+100)*(s+13)); % L(s) for 5.13

L = L1; % Choose root locus

%% Plot root locus
figure
rlocus(L)

%% Find zeta = 0.707

% See math in homework, coefficients of characteristic equation such that s = sigma + sigma*j (zeta = sqrt(2)/2)
% Characteristic equation: s^5 + 13s^4 + (100+K)s^3 + (1300+K)s^2 + 81Ks + 81K = 0
polynom = [16 112 104 -912 32400 226800 210600]; 

sigRaw = roots(polynom);
sigmas = sigRaw(imag(sigRaw) == 0) % Only interested in real sigmas

K = sigmas.^3.*(4*sigmas.^2+52*sigmas+200)./(-2*sigmas.^3+81*sigmas+81)

%% Simulate step response with these K's

for k = 1:length(K)
    T(k) = feedback(K(k).*L,1); % Unity feedback of K*L(s)

    figure
    step(T(k))
    titleText = sprintf("Step response for K = %.5f", K(k));
    title(titleText);
end

T