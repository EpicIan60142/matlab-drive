%% ASEN 3200 HW O-5 Problem 7.12 Script
%   By: Ian Faber, 11/28/2022

%% Housekeeping
clc; clear; close all

%% Constants and givens
mu = 398600; % km^3/s^2, Earth
rStation = 6600; % km
TStation = 2*pi*sqrt(rStation^3/mu); % sec
T = TStation/3; % sec

r0 = [1; 1; 1]; % km
v0 = [0; 0; 0.005]; % km/s
vf = [0; 0; 0]; % m/s

n = (2*pi)/TStation ; % rad/s

size = 35; % Marker size for 3D plot

%% Relative Motion matrices
Phi_rr = @(t,n) [
                    4 - 3*cos(n*t),         0,      0;
                    6*(sin(n*t) - n*t),     1,      0;
                    0,                      0,      cos(n*t)
                ];

Phi_rv = @(t,n) [
                    (1/n)*sin(n*t),         (2/n)*(1-cos(n*t)),     0;
                    (2/n)*(cos(n*t) - 1),   (4/n)*sin(n*t) - 3*t,   0;
                    0,                      0,                      (1/n)*sin(n*t)
                ];

Phi_vr = @(t,n) [
                    3*n*sin(n*t),           0,      0;
                    6*(n*cos(n*t)-n),       0,      0;
                    0,                      0,      -n*sin(n*t)
                ];

Phi_vv = @(t,n) [
                    cos(n*t),       2*sin(n*t),         0;
                    -2*sin(n*t),    4*cos(n*t) - 3,     0;
                    0,              0,                  cos(n*t)
                ];


%% Problem

v0plus = -(Phi_rv(T,n)^-1)*Phi_rr(T,n)*r0; % km/s
deltaV1 = norm(v0plus - v0) % m/s

% n = (2*pi)/T;

vfminus = Phi_vr(T,n)*r0 + Phi_vv(T,n)*v0plus; % km/s
deltaV2 = norm(vf - vfminus) % m/s

deltaVTotal = deltaV1 + deltaV2 % m/s

%% Plot Trajectory
stop = TStation; % Full sim: TStation, Up to rendezvous: T

t = (0:1:stop);

for k = 1:length(t)
    r(k,:) = Phi_rr(t(k),n)*r0 + Phi_rv(t(k),n)*v0plus;
end

figure
hold on
grid on
grid minor
title("Problem 7.12 Chaser Vehicle Trajectory")
plot3(r(:,1), r(:,2), r(:,3))
scatter3(r0(1), r0(2), r0(3), size, "green", 'filled')
scatter3(0, 0, 0, size, "black", 'filled')
scatter3(r(end,1), r(end,2), r(end,3), size, "red", "filled")
xlabel("Radial (km)")
ylabel("In-track (km)")
zlabel("Cross-track (km)")
view([30 35])

legend("Trajectory", "Start", "Rendezvous", "End")



