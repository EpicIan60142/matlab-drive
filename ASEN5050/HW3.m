%% ASEN 5050 HW 3 Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

addpath("utilities\")

%% Problem 1a
% Setup
mu = 3.794e7; % km^3/s^2
rSaturn = 60268; % km

% Initial time t1
r1 = [-720000; 670000; 310000]; % km
r = norm(r1);
v1 = [2.160; -3.360; 0.620]; % km/s
v = norm(v1);

% Calculate useful vectors and quantities
eVec = ((v^2 - (mu/r))*r1 - dot(r1,v1)*v1)/mu
e = norm(eVec)

a = (-mu)/(v^2 - ((2*mu)/r))

hVec = cross(r1, v1)
h = norm(hVec)

% Solve for cos(trueAnom) at impact
cosTA = (1/e)*(((a*(1-e^2))/rSaturn) - 1)
sinTA = 1-cosTA^2

% Solve for vectors at impact in r theta h frame
r2 = [rSaturn; 0; 0]
v2 = [(mu/h)*e*sinTA; (mu/h)*(1+e*cosTA); 0]

% Get orbital elements and find DCM
TA = acos(cosTA)*sign(dot(r2, v2))

nHat = cross([0; 0; 1], hVec/h)
RAAN = acos(dot(nHat, [1;0;0]))*sign(dot(nHat,[0;1;0]))

inc = acos(dot(hVec/h, [0;0;1]))

argPeri = acos(dot(nHat, eVec/e))*sign(dot(eVec/e,[0;0;1]))

theta = TA + argPeri

C = EA2DCM([-theta, -inc, -RAAN], [3, 1, 3])

% Calculate position and velocity vectors in inertial frame
rXYZ = C*r2
vXYZ = C*v2

%% Problem 2c
% Setup
mu = 42828.4; % km^3/s^2

e = 0.45454;
a = 6463.8; % km
inc = 74.924; % deg
RAAN = 1.241; % deg
argPeri = 353.31; % deg
TA = 199.38; % deg

% Calculate theta
theta = mod(TA + argPeri,360); % deg

% Calculate h, r, vr, and vtheta
h = sqrt(mu*a*(1-e^2));
r = (a*(1-e^2))/(1+e*cosd(TA));
vr = (mu/h)*e*sind(TA); 
vtheta = (mu/h)*(1+e*cosd(TA));

r1 = [r; 0; 0]
v1 = [vr; vtheta; 0]

C = EA2DCM(deg2rad([-theta, -inc, -RAAN]), [3,1,3])

rVec = C*r1
vVec = C*v1


