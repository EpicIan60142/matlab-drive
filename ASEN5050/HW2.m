%% ASEN 5050 HW 2 Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all
format long
addpath("utilities\")

%% let's do this bullshit
mu = 4.305e4;
r_vec = [3.62067e3; -3.19925e2; -4.20645e2];
r = norm(r_vec);
v_vec = [-4.28843e-1; -3.00176e-2; -3.39801];
v = norm(v_vec);

e_vec = ((v^2 - mu/r)*r_vec - dot(r_vec, v_vec)*v_vec)/mu;
e = norm(e_vec);

h_vec = cross(r_vec, v_vec);
h = norm(h_vec);

n_vec = cross([0; 0; 1], h_vec);
n = norm(n_vec);

inc = acosd(dot(h_vec, [0; 0; 1])/h)
RAAN = acosd(dot(n_vec,[1; 0; 0])/n)*sign(dot(n_vec, [0; 1; 0]))
argPeri = acosd(dot(n_vec, e_vec)/(n*e))*sign(dot(e_vec,[0; 0; 1]))
trueAnom = acosd(dot(r_vec, e_vec)/(r*e))*sign(dot(r_vec, v_vec))

theta = argPeri + trueAnom

C = EA2DCM(deg2rad([-theta, -inc, -RAAN]), [3,1,3])

r_vec_2 = C'*r_vec

v_vec_2 = C'*v_vec

% Now, e_hat = r_hat
e_hat = e_vec/e
h_hat = h_vec/h
t_hat = cross(h_hat, e_hat)

r_vec_3 = r*e_hat

vTheta = (mu/h)*(1+e);
v_vec_3 = vTheta*t_hat

format short