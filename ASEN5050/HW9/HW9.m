%% ASEN 5050 HW 9 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
muEarth = 3.986004415e5; % km^3/s^2
muMars = 4.305e4; % km^3/s^2

rEarth = 6378.1363; % km
rMars = 3397.2; % km

Prot_Mars = 1.02595675; % days

r0 = 10000; % km

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

%% Problem 1
fprintf("--- Problem 1 ---\n")
x0 = 5; % m
xDot0 = 0; % m/s
z0 = 0 % m -> Same orbit plane
zDot0 = 0 % m/s -> Same orbit plane
yMax = 15; % m

n = sqrt(muEarth/(r0^3));

yDot0 = -2*n*x0 % -> Drift terms suppressed
y0 = yMax - 2*x0

%% Problem 2
fprintf("--- Problem 2 ---\n")
format shortE
r0_vec = [2;2;0]; % m
v0minus = [-0.03; 0.01; 0.05]; % m/s
r1_vec = [-2;2;0]; % m

T = 2*pi/n; % sec
t1 = T/2;

v0plus = (Phi_rv(t1,n)^-1)*(r1_vec - Phi_rr(t1,n)*r0_vec)
dv1 = norm(v0plus - v0minus)

v1minus = Phi_vr(t1,n)*r0_vec + Phi_vv(t1,n)*v0plus
dv2 = norm(zeros(size(v1minus)) - v1minus)

    % Plot Trajectory
t = 0:t1;

markerSize = 25;

for k = 1:length(t)
    r(:,k) = Phi_rr(t(k),n)*r0_vec + Phi_rv(t(k),n)*v0plus;
    v(:,k) = Phi_vr(t(k),n)*r0_vec + Phi_vv(t(k),n)*v0plus;
end

figure
hold on; grid on; grid minor;
title("Problem 2 CubeSat Trajectory")
traj = plot3(r(1,:), r(2,:), r(3,:));
plot3(r0_vec(1), r0_vec(2), r0_vec(3), 'g.', 'MarkerSize', markerSize)
plot3(0, 0, 0, 'k.', 'MarkerSize', markerSize)
plot3(r(1,end), r(2,end), r(3,end), 'r.', 'MarkerSize', markerSize)
xlabel("Radial [m]"); ylabel("Along-track [m]"); zlabel("Cross-track [m]")
xlim([-5 2.5]); ylim([-2.5 5]); zlim([-2 2]); view([-90 90])

legend("Trajectory", "Start", "Primary s/c", "End")

labels = ["Xdot", "Ydot", "Zdot"];
values = [round(v(1,:),7); round(v(2,:),7); round(v(3,:),7)];
for k = 1:length(labels)
    extraRows(k) = dataTipTextRow(labels(k), values(k,:));
    traj.DataTipTemplate.DataTipRows(end+length(extraRows)) = extraRows(k);
end

points = [1, size(r,2)];
for k = 1:length(points)
    datatip(traj, r(1,points(k)), r(2,points(k)), r(3,points(k)), 'location', 'southwest');
end

%% Problem 3
fprintf("--- Problem 3 ---\n")
    % a
omegaMars = (2*pi)/(Prot_Mars*86400); % rad/s

dLong = deg2rad(-210);
a = (muMars*(dLong/(-2*pi*omegaMars))^2)^(1/3)

format default

