%% ASEN 5050 HW 7 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
muSun = 1.32712428e11; % km^3/s^2
muEarth = 3.986004415e5; % km^3/s^2
muSaturn = 3.794e7; % km^3/s^2
muJupiter = 1.268e8; % km^3/s^2

aEarth = 1.0000010178; % AU
aSaturn = 9.554909595; % AU
aJupiter = 5.202603191; % AU
AU2km = 149597870.7; % 1 AU = 149,597,870.7 km

rSaturn = 60268; % km
rJupiter = 71492; % km

G = 6.673e-20; % km^3/kg-s^2

%% Problem 1 setup
r_a = 1.8e6; % km
r_p = 600000; % km

mTitan = 1.3455e23; % kg
aTitan = 1.221830e6; % km
rTitan = 2575; % km

%% Problem 1a
e_in = (r_a - r_p)/(r_a + r_p);
a_in = 0.5*(r_a + r_p);

TAin_deg = -acosd((1/e_in)*((a_in*(1-e_in^2)/aTitan) - 1)) % Negative since flying from apoapsis to periapsis

h = sqrt(muSaturn*a_in*(1-e_in^2));

vr = (muSaturn/h)*e_in*sind(TAin_deg);
vtheta = (muSaturn/h)*(1+e_in*cosd(TAin_deg));

Vin = [vr; vtheta] % rhat; thetahat

VTitan = [0; sqrt(muSaturn/aTitan)];

VinfIn = Vin - VTitan

%% Problem 1b
r_ph = 3000; % km

muTitan = G*mTitan;

vinfIn = norm(VinfIn);

a_h = -muTitan/(vinfIn^2);
e_h = 1-(r_ph/a_h);

delta_deg = 2*asind(1/e_h)

%% Problem 1e
vin = norm(Vin);
vTitan = norm(VTitan);

beta_1_deg = acosd((vin^2 - vinfIn^2 - vTitan^2)/(-2*vinfIn*vTitan));
beta_2_deg = beta_1_deg - delta_deg;

vinfOut = vinfIn;
vout = sqrt(vinfOut^2 + vTitan^2 - 2*vinfOut*vTitan*cosd(beta_2_deg));

fpaOut_deg = -acosd((vinfOut^2 - vTitan^2 - vout^2)/(-2*vTitan*vout)); % Negative since flying from apoapsis to periapsis

Vout = [vout*sind(fpaOut_deg); vout*cosd(fpaOut_deg)]

%% Problem 1f
rout = aTitan;
Rout = [aTitan; 0];

a_out = -muSaturn/(vout^2 - ((2*muSaturn)/rout))

eVec_out = ((vout^2 - muSaturn/rout)*Rout - dot(Rout, Vout)*Vout)/muSaturn
e_out = norm(eVec_out)

TAout_deg = -acosd((1/e_out)*((a_out*(1-e_out^2)/rout) - 1)) % Negative since flying from apoapsis to periapsis

r_aOut = a_out*(1+e_out);
r_pOut = a_out*(1-e_out);

%% Problem 1g
dV_eq = Vout - Vin
dv_eq = norm(dV_eq)


