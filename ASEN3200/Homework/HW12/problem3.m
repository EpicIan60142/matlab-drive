%% Housekeeping
clc; clear; close all

%% Constants
mu = 132712440018; % km^3/s^2, Sun

r1 = [-1.054e8; 1.579e8; -1.520e5]; % km
d1 = norm(r1);
r2 = [-1.461e8; 1.081e8; -2.265e5]; % km
d2 = norm(r2);
r3 = [-1.652e8; 4.254e7; -2.673e5]; % km
d3 = norm(r3);

%% Problem

% 3.a
N = d1*cross(r2,r3) + d2*cross(r3,r1) + d3*cross(r1,r2);
n = norm(N);
D = cross(r1,r2) + cross(r2,r3) + cross(r3,r1);
d = norm(D);
S = r1*(d2-d3) + r2*(d3-d1) + r3*(d1-d2);

v2 = sqrt(mu/(n*d))*((cross(D,r2)/d2)+S) % km/s

% 3.b
h = cross(r2, v2);
inc = rad2deg(acos(h(3)/norm(h)))

a = -mu/(norm(v2)^2 - 2*(mu/d2))

ecc = (cross(v2, h)/mu) - (r2/d2);
e = norm(ecc)

% 3.c
theta = rad2deg(acos(((a*(1-e^2))/(d2*e)) - (1/e)))
