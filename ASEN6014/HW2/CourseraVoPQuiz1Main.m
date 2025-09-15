%% ASEN 6014 Lagrange Matrix Coursera Quiz Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
const.k1 = 3;
const.k3 = 1;

e0 = [1; 0];

tspan = 0:0.01:10;

%% Integrate

[t, X] = ode45(@(t,X)springInvariantEOM(t,X,const), tspan, e0, odeset('RelTol',1e-12,'AbsTol',1e-12));

eEnd = X(end,:)'

e1 = X(:,1);
e2 = X(:,2);

x = e1.*cos(const.k1*t) + (e2/const.k1).*sin(const.k1*t);

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
title(tl, "Evolution of x")
nexttile;
    hold on; grid on;
    plot(t, x, 'b');
    xlabel("Time [sec]"); ylabel("x")

