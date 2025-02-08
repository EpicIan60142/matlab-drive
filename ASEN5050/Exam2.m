%% ASEN 5050 Exam 2 Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
muSun = 1.32712428e11; % km^3/s^2
muMoon = 4.902799e3; % km^3/s^2
muJupiter = 1.268e8; % km^3/s^2

aJupiter = 5.202603191; % AU

rMoon = 1738; % km
rJupiter = 71492; % km

AU2km = 149597870.7; % 1 AU = 149,597,870.7 km

g0 = 9.81; % m/s^2

%% Problem 1a
fprintf("--- Problem 1a ---\n")
rp_i = 1850; % km
ra_i = 5400; % km
r_i = 3500; % km

rp_f = 3000; % km
ra_f = 5400; % km

a_i = 0.5*(ra_i + rp_i)
e_i = (ra_i - rp_i)/(ra_i + rp_i)
h_i = sqrt(muMoon*a_i*(1-e_i^2))
TAi_deg = -acosd((1/e_i)*(((a_i*(1-e_i^2))/r_i) - 1)) % Negative since moving from apoapsis to periapsis
V_i = [(muMoon/h_i)*e_i*sind(TAi_deg); (muMoon/h_i)*(1+e_i*cosd(TAi_deg))]

a_f = 0.5*(ra_f + rp_f)
e_f = (ra_f - rp_f)/(ra_f + rp_f)
h_f = sqrt(muMoon*a_f*(1-e_f^2))
r_f = r_i; % Impulsive maneuver happens instantly
TAf_deg = -acosd((1/e_f)*(((a_f*(1-e_f^2))/r_f) - 1)) % Can be positive or negative - negative gives smallest dV magnitude
V_f = [(muMoon/h_f)*e_f*sind(TAf_deg); (muMoon/h_f)*(1+e_f*cosd(TAf_deg))]

dV = V_f - V_i
dv = norm(dV)

%% Problem 1b
fprintf("--- Problem 1b ---\n")
dArgPeri = TAi_deg - TAf_deg

%% Problem 2a
fprintf("--- Problem 2a ---\n")
muX = 4000; % km^3/s^2

a_in = 1.8e6; % km
V_in = [2.2789; 5.8841]; % km/s
v_in = norm(V_in)

r_in = muJupiter/((v_in^2)/2 + muJupiter/(2*a_in))

%% Problem 2b
fprintf("--- Problem 2b ---\n")
rX = r_in;

VX = [0; sqrt(muJupiter/rX)]

Vinf_in = V_in - VX

%% Problem 2c
fprintf("--- Problem 2c ---\n")
v_out = 7.6759; % km/s
phi_out_deg = 20.91; % deg

V_out = [v_out*sind(phi_out_deg); v_out*cosd(phi_out_deg)]

%% Problem 2d
fprintf("--- Problem 2d ---\n")

Vinf_out = V_out - VX

%% Problem 2h
fprintf("--- Problem 2h ---\n")
vinf_in = norm(Vinf_in);
vinf_out = norm(Vinf_out);

a_h = -muX/(vinf_in^2);

delta_deg = acosd(dot(Vinf_in, Vinf_out)/(vinf_in*vinf_out));

e_h = 1/sind(delta_deg/2);

rp_h = a_h*(1-e_h)

%% Problem 2i
fprintf("--- Problem 2i ---\n")

dV_eq = V_out - V_in;
dv_eq = norm(dV_eq)

%% Problem 2j
fprintf("--- Problem 2j ---\n")
m_i = 2500; % kg
m_p = 500; % kg
Isp = 300; % sec

dv = 1000*dv_eq; % km/s -> m/s
m_p_req = m_i*(1-exp(-dv/(Isp*g0)))

%% Problem 3
fprintf("--- Problem 3 ---\n")
Vin = [3.2476; 13.1175]; % km/s
Vout = [3.3426; 13.7356]; % km/s

VJ = [0; sqrt(muSun/(AU2km*aJupiter))];

Vinf_in = Vin - VJ;
vinf_in = norm(Vinf_in)

Vinf_out = Vout - VJ;
vinf_out = norm(Vinf_out)
