%% ASEN 5050 HW 5 Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Constants
muMars = 4.305e4; % km^3/s^2 
muMoon = 4902.799; % km^3/s^2
muSun = 1.32712428e11; % km^3/s^2
muSaturn = 3.794e7; % km^3/s^2

RMars = 3397.2; % km
RMoon = 1738; % km

aEarth = 1.0000010178; % AU
aSaturn = 9.554909595; % AU

AU2km = 149597870.7; % 1 AU = 149,597,870.7 km

%% Problem 1
a1 = 7045; % km
e1 = 0.23;
TA1_deg = -142; 

% Part a
fprintf("--- Problem 1 part a ---")
h1 = sqrt(muMoon*a1*(1-e1^2));
vr1 = (muMoon/h1)*e1*sind(TA1_deg);
vtheta1 = (muMoon/h1)*(1+e1*cosd(TA1_deg));

r1 = (a1*(1-e1^2))/(1+e1*cosd(TA1_deg));

v1 = [vr1; vtheta1]

% Part c
fprintf("--- Problem 1 part c ---")
dv = [0.3; -0.1];
v2 = v1 + dv

% Part d
fprintf("--- Problem 1 part d ---")
v = norm(v2);
fpa2 = asind(v2(1)/v);
r2 = r1;
a2 = -muMoon/(v^2 - (2*muMoon)/r2)
specEng2 = -muMoon/(2*a2);
h2 = r2*v2(2);

e2 = sqrt(1+2*((h2^2*specEng2)/(muMoon^2)))

TA2_deg = acosd((1/e2)*(((h2*v2(2))/muMoon) - 1))*sign(fpa2)

% Part e
fprintf("--- Problem 1 part e ---")
dArgPeri = TA1_deg - TA2_deg

% Part f
fprintf("--- Problem 1 part f ---")
mi = 1224; % kg
g0 = 9.81; % m/s^2
Isp = 212; % sec
magDV = 1000*norm(dv); % km/s -> m/s

mp = mi*(1-exp((-magDV)/(g0*Isp)))

%% Problem 2
r1 = 6500; % km
E1 = pi/2; % rad
rp1 = 5915; % km
rp2 = 5712; % km
ra2 = 7888; % km

% Part a
fprintf("--- Problem 2 part a ---")
a1 = r1/(1-e1*cos(E1));
e1 = 1-(rp1/a1);
h1 = sqrt(muMars*a1*(1-e1^2));

TA1_deg = 2*atand(sqrt((1+e1)/(1-e1))*tan(E1/2));

vr1 = (muMars/h1)*e1*sind(TA1_deg);
vtheta1 = (muMars/h1)*(1+e1*cosd(TA1_deg));

v1 = [vr1; vtheta1]

% Part b
fprintf("--- Problem 2 part b ---")
a2 = 0.5*(ra2 + rp2);
e2 = (ra2 - rp2)/(ra2 + rp2);
h2 = sqrt(muMars*a2*(1-e2^2));
r2 = r1;
TA2_deg = acosd((1/e2)*(((a2*(1-e2^2))/r2)-1));

vr2 = (muMars/h2)*e2*sind(TA2_deg);
vtheta2 = (muMars/h2)*(1+e2*cosd(TA2_deg));

v2 = [vr2; vtheta2];

dv = v2 - v1;
magDV = norm(dv)

%% Problem 3
ai = aEarth*AU2km; % AU -> km
af = aSaturn*AU2km; % AU -> km

% Part a
fprintf("--- Problem 3 part a ---")
at_HT = 0.5*(ai + af); % km

v1i_HT = sqrt(muSun/ai);
v1f_HT = sqrt(((2*muSun)/ai) - (muSun/at_HT));
v2i_HT = sqrt(((2*muSun)/af) - (muSun/at_HT));
v2f_HT = sqrt(muSun/af);

dv1_HT = v1f_HT - v1i_HT;
dv2_HT = v2f_HT - v2i_HT;

dvTotal_HT = abs(dv1_HT) + abs(dv2_HT)

TOF_HT = pi*sqrt((at_HT^3)/muSun);
TOF_HT_yr = TOF_HT/31536000 % sec -> years

% Part b
fprintf("--- Problem 3 part b ---")
alpha_HT = sqrt(muSun/(af^3))*TOF_HT
phi_HT = pi - alpha_HT;
phi_HT_deg = rad2deg(phi_HT)

% Part c
fprintf("--- Problem 3 part c ---")
rB = 11*AU2km; % AU -> km
at1_BE = 0.5*(ai + rB);
at2_BE = 0.5*(af + rB);

v1i_BE = sqrt(muSun/ai);
v1f_BE = sqrt(((2*muSun)/ai) - (muSun/at1_BE));
v2i_BE = sqrt(((2*muSun)/rB) - (muSun/at1_BE));
v2f_BE = sqrt(((2*muSun)/rB) - (muSun/at2_BE));
v3i_BE = sqrt(((2*muSun)/af) - (muSun/at2_BE));
v3f_BE = sqrt(muSun/af);

dv1_BE = v1f_BE - v1i_BE;
dv2_BE = v2f_BE - v2i_BE;
dv3_BE = v3f_BE - v3i_BE;

dvTotal_BE = abs(dv1_BE) + abs(dv2_BE) + abs(dv3_BE)

TOF_BE = pi*(sqrt((at1_BE^3)/muSun) + (sqrt((at2_BE^3)/muSun)));
TOF_BE_yr = TOF_BE/31536000

% Part d
fprintf("--- Problem 3 part d ---")
k1 = rB/ai
k2 = af/ai

