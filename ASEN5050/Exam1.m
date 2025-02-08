%% ASEN 5050 Exam 1 Script
% By: Ian Faber
%   - MATLAB was allowed on this exam for working!


%% Housekeeping
clc; clear; close all

%% Setup
mu = 4.305e4; % km^3/s^2
rMars = 3397.2; % km

%% Problem 2a
R3 = [-7.665e3; 6.5468e3; -4.574e2];
r3 = norm(R3)
V3 = [1.6334; 0.1226; -1.9455];
v3 = norm(V3)

specEng = (v3^2)/2 - (mu/r3)
a = -mu/(2*specEng)

e_vec = ((v3^2 - mu/r3)*R3 - dot(R3, V3)*V3)/mu
e = norm(e_vec)

h_vec = cross(R3, V3)
h = norm(h_vec)

n_vec = cross([0; 0; 1], h_vec)
n = norm(n_vec);

inc = acosd(dot(h_vec, [0;0;1])/h)
RAAN = acosd(dot(n_vec,[1; 0; 0])/n)*sign(dot(n_vec, [0; 1; 0]))
argPeri = acosd(dot(n_vec, e_vec)/(n*e))*sign(dot(e_vec,[0; 0; 1]))
TA3_deg = acosd(dot(R3, e_vec)/(r3*e))*sign(dot(R3, V3))

%% Problem 2b
TA3 = deg2rad(TA3_deg)
E3 = 2*atan(sqrt((1-e)/(1+e))*tan(TA3/2))

T = (2*pi*sqrt((a^3)/mu))

t3 = (T/(2*pi))*(E3 - e*sin(E3))
dt = 2*3600; % 2 hours
t4 = (t3 + dt)

E4 = solveKeplersEq(t4, a, e, mu)

TA4 = 2*atan(sqrt((1+e)/(1-e))*tan(E4/2))
TA4_deg = rad2deg(TA4)

dTA = TA4 - TA3
dTA_deg = rad2deg(dTA)

p = a*(1-e^2)
r4 = p/(1+e*cos(TA4))

f = 1 - (r4/p)*(1 - cos(dTA))
g = ((r4*r3)/sqrt(mu*p))*sin(dTA)
fDot = sqrt(mu/p)*tan(dTA/2)*(((1-cos(dTA))/p) - (1/r4) - (1/r3))
gDot = 1 - (r3/p)*(1 - cos(dTA))

R4 = f*R3 + g*V3
V4 = fDot*R3 + gDot*V3

