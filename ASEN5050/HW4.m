%% ASEN 5050 HW 4 Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

addpath("utilities\")

%% Problem 1
% Setup
fprintf("------ Problem 1 ------")
mu = 3.794e7; % km^3/s^2

r1 = [-720000; 670000; 310000]; % km
r0 = norm(r1);
v1 = [2.16; -3.36; 0.62]; % km/s
v0 = norm(v1);

% Calculate helper quantities
h = norm(cross(r1, v1));
p = (h^2)/mu;
eng = (v0^2)/2 - (mu/r0);
e = sqrt(1 + ((2*(h^2)*eng)/(mu^2)));

TA1 = acosd((1/e)*((p/r0) - 1))*sign(dot(r1, v1));

TA2 = -13.17; % deg
r2 = p/(1+e*cosd(TA2));

dTA = TA2 - TA1;

% Find f and g function values
f = 1 - (r2/p)*(1 - cosd(dTA));
g = ((r2*r0)/sqrt((mu*p)))*sind(dTA);
fDot = sqrt(mu/p)*tand(dTA/2)*(((1 - cosd(dTA))/p) - (1/r2) - (1/r0));
gDot = 1 - (r0/p)*(1 - cosd(dTA));

% Calculate r2 and v2 (1a)
r2 = f*r1 + g*v1
v2 = fDot*r1 + gDot*v1

% Get helper quantities for 1b
a = p/(1-e^2);

T = 2*pi*sqrt((a^3)/mu);

% Solve for eccentric anomalies
E2 = 2*atan(sqrt((1-e)/(1+e))*tand(TA2/2));
E1 = 2*atan(sqrt((1-e)/(1+e))*tand(TA1/2));

% Calculate time elapsed (1b)
t12 = (T/(2*pi))*((E2 - e*sin(E2)) - (E1 - e*sin(E1)));
t12 = t12/3600 % convert from sec to hr

%% Problem 2
fprintf("------ Problem 2 ------")
% Setup
mu = 1.268e8; % km^3/s^2

R1 = [5.35295e6; 7.053778e5; -4.0597e5]; % km
r1 = norm(R1);
V1 = [-4.164248; 1.96369; 3.191257e-1]; % km/s
v1 = norm(V1);

% Calculate helper quantities
e_vec = ((v1^2 - (mu/r1))*R1 - dot(R1, V1)*V1)/mu;
e = norm(e_vec);

a = -mu/(v1^2 - ((2*mu)/r1));

% Calculate true and eccentric anomaly at t1 (2a)
TA1 = acos((1/e)*((a*(1-e^2)/r1) - 1))*sign(dot(R1, V1));
TA1_deg = rad2deg(TA1)

E1 = 2*atan(sqrt((1-e)/(1+e))*tan(TA1/2))

% Calculate helper quantities for 2b
h = cross(R1, V1);

nHat = cross([0;0;1], h)/norm(cross([0;0;1], h));

inc = rad2deg(acos(dot(h/norm(h),[0;0;1])));
RAAN = rad2deg(acos(dot(nHat,[1;0;0]))*sign(dot(nHat, [0;1;0])));
argPeri = rad2deg(acos(dot(nHat,e_vec/e))*sign(dot(e_vec, [0;0;1])));

% Calculate true and eccentric anomaly at t2 (2b)
TA2_deg = 180 - argPeri
TA2 = deg2rad(TA2_deg);

E2 = 2*atan(sqrt((1-e)/(1+e))*tan(TA2/2))

% Calculate period for 2c
T = 2*pi*sqrt((a^3)/mu);
T_days = T/(60*60*24);

% Calculate time elapsed from t1 to t2 (2c)
t12 = (T_days/(2*pi))*((E2 - e*sin(E2)) - (E1 - e*sin(E1)))

% Calculate helper quantities for 2d
p = a*(1-e^2);
r2 = p/(1+(e*cos(TA2)));
dTA = TA2 - TA1;
dTA_deg = rad2deg(dTA);

% Calculate f and g function values at t2
f = 1 - (r2/p)*(1 - cos(dTA));
g = ((r2*r1)/sqrt(mu*p))*sin(dTA);
fDot = sqrt(mu/p)*tan(dTA/2)*(((1-cos(dTA))/p) - (1/r2) - (1/r1));
gDot = 1 - (r1/p)*(1 - cos(dTA));

% Calculate position and velocity vectors at t2 (2d)
R2 = f*R1 + g*V1
V2 = fDot*R1 + gDot*V1

check = -r2*nHat;

% Setup for 2f
dt = 20*(24*60*60); % Time in seconds since t1

t1 = sqrt((a^3)/mu)*(E1 - e*sin(E1));
t3 = t1 + dt;

% Solve for E3 and TA3 at t3 (2f)
E3 = solveKeplersEq(t3, a, e, mu)
TA3 = 2*atan(sqrt((1+e)/(1-e))*tan(E3/2));
TA3_deg = rad2deg(TA3)


