%% ASEN 5010 HW 3 Concept Check 20 Script
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath('..\Utilities\')

%% Problem
sig0 = [0.4; 0.2; -0.1];

t = (0:0.001:60)';

[time, sig] = ode45(@(t, state)MRPEOM(t, state), t, sig0);

% sigNorms = zeros(length(sig),1);
for k = 1:length(sig)
    sigNorm = norm(sig(k,:));

    if sigNorm > 1
        sig(k,:) = -sig(k,:)/sigNorm^2;
    end
end

timeInstant = 42;

answer = norm(sig(time == timeInstant,:))

% b0 = zeros(length(sig),1);
% 
% for k = 1:length(sig)
%     b0(k) = 1./sqrt(1 + sig(k,:)*sig(k,:)');
%     phi(k) = 2*acos(b0(k));
% end

figure
hold on
sgtitle("Modified Rodriquez Parameter Evolution")

subplot(3,1,1)
plot(time, sig(:,1));
xlabel("time [sec]")
ylabel("\sigma_1")

subplot(3,1,2)
plot(time, sig(:,2));
xlabel("time [sec]")
ylabel("\sigma_2")

subplot(3,1,3)
plot(time, sig(:,3));
xlabel("time [sec]")
ylabel("\sigma_3")