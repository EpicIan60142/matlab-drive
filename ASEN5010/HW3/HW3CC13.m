%% ASEN 5010 HW 3 Concept Check 13 Script
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath('..\Utilities\')

%% Problem
q0 = [0.4; 0.2; -0.1];

t = (0:0.001:60)';

[time, q] = ode45(@(t, state)CRPEOM(t, state), t, q0);

timeInstant = 42;

answer = norm(q(time == timeInstant,:))

b0 = zeros(length(q),1);

for k = 1:length(q)
    b0(k) = 1./sqrt(1 + q(k,:)*q(k,:)');
    phi(k) = 2*acos(b0(k));
end

figure
hold on
sgtitle("Classical Rodriquez Parameter Evolution")

subplot(3,1,1)
plot(time, q(:,1));
xlabel("time [sec]")
ylabel("q_1")

subplot(3,1,2)
plot(time, q(:,2));
xlabel("time [sec]")
ylabel("q_2")

subplot(3,1,3)
plot(time, q(:,3));
xlabel("time [sec]")
ylabel("q_3")
