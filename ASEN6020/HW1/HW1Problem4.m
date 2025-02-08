%% ASEN 6020 HW 1 Problem 4
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup cost

J_3B = @(r, i1, i3) sqrt(1 + ((2*r)./(1+r)) - 2*sqrt(2*r./(1+r)).*cos(i1)) + 2*sqrt(2/(r.*(1+r))).*sin(((pi/2)-i1-i3)/2) + sqrt(1 + ((2*r)./(1+r)) - 2*sqrt(2*r./(1+r)).*cos(i3));

%% Create contour plot
r = 10;

i1 = deg2rad(0:0.5:90);
i3 = deg2rad(0:0.5:90);

[I1_full, I3_full] = meshgrid(i1, i3);

idx = pi/2 - I1_full - I3_full <= 0;
I1_full(idx) = NaN;
I1 = I1_full;
I3_full(idx) = NaN;
I3 = I3_full;


J = J_3B(r, I1, I3);

[minCost, minIdx] = min(J,[],'all');
i1_min = I1(minIdx);
i3_min = I3(minIdx);

figure; hold on; grid on;
title("Cost of Dog-Leg Bi-elliptic Plane Change Maneuver")
levels = logspace(log10(minCost), log10(2.5), 300);
contour(rad2deg(I1), rad2deg(I3), J, levels); c = colorbar;
ax = plot3(rad2deg(i1_min), rad2deg(i3_min), minCost, 'k.', 'MarkerSize', 10);
c.Label.String = "Cost"; datatip(ax, rad2deg(i1_min), rad2deg(i3_min), 'Location', 'northeast');
xlabel("\Deltai_1 [deg]"); ylabel("\Deltai_3 [deg]");
