%% ASEN 6020 HW 1 Problem 5 Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
costFun = @(eta, l, i) sqrt((2*l)/(1+l) + 1 - 2*sqrt((2*l)/(1+l))*cos(eta*i)) + sqrt((4/(l*(1+l)))*(1 - cos((1-eta)*i)));
fun = @(eta, l, i) (sqrt((2*l)/(1+l))*i*sin(eta*i))/sqrt(((2*l)/(1+l)) + 1 - 2*sqrt((2*l)/(1+l))*cos(eta*i)) - ...
                   ((2/(l*(l+1)))*i*sin((1-eta)*i))/sqrt(((4*l)/(l+1))*(1 - cos((1-eta)*i)));

% f = @(x) fun(x, l, i); % x = [eta; l; i]

%% Solve for etaStar given an l and delta_i
    % Set test params
l_nom = 10;
i_nom = deg2rad(90);

    % Test function
etaStar = fsolve(@(eta)fun(eta, l_nom, i_nom), 0)

%% Investigate how etaStar changes as a function of delta_i and l
    % Disable fsolve summary
opt = optimoptions("fsolve","Display","none");

    % Set test regime
l_test = 1:100;
i_test = deg2rad(4:2:90);

    % Calculate etaStar for each l and delta_i
etaStar = [];
cost = [];
for k = 1:length(i_test)
    i_nom = i_test(k);
    eta_part = [];
    cost_part = [];
    for kk = 1:length(l_test)
        l_nom = l_test(kk);
        eta = fsolve(@(eta)fun(eta, l_nom, i_nom), 0.01, opt);
        eta_part = [eta_part; eta];
        cost_part = [cost_part; costFun(eta, l_nom, i_nom)];
    end
    etaStar = [etaStar, eta_part];
    cost = [cost, cost_part];
end

    % Plot results
figure
grid on; hold on;
title("\eta* vs. l for varying \Deltai")
x = l_test.*ones(length(i_test), length(l_test));
y = etaStar;
z = rad2deg(i_test).*ones(length(l_test), length(i_test));
contour(x', y, z, rad2deg(i_test)); c = colorbar;
c.Label.String = "\Deltai [deg]"; colormap("turbo")
xlabel("l"); ylabel("\eta*");
set(gca, 'yscale', 'log')









