%% ASEN 6014 Lagrangian Brackets Coursera Quiz Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
const.m = 1;
const.g = 9.81;
const.k = 3;

w = sqrt(const.k/const.m);

e0 = [1; 0];

tspan = 0:0.01:10;

%% Integrate

[t, X] = ode45(@(t,X)gravityInvariantEOM(t,X,const), tspan, e0, odeset('RelTol',1e-12,'AbsTol',1e-12));

eEnd = X(end,:)'

e1 = X(:,1);
e2 = X(:,2);

y = e1.*cos(w*t) + (e2/w).*sin(w*t);

%% Plot
figure; tl = tiledlayout(2,1);
title(tl, "Evolution of e")
nexttile;
    hold on; grid on;
    plot(t, X(:,1), 'b');
    xlabel("Time [sec]"); ylabel("e_1")
nexttile;
    hold on; grid on;
    plot(t, X(:,2), 'r');
    xlabel("Time [sec]"); ylabel("e_2")

figure; tl = tiledlayout(2,1);
title(tl, "Evolution of y")
nexttile;
    hold on; grid on;
    plot(t, y, 'b');
    xlabel("Time [sec]"); ylabel("y")

