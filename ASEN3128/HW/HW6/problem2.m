%% Housekeeping
clc; clear; close all;

%% Constants
% Matrices
A = [-0.8771 0.1146 0; -0.2797 6.976e-3 -1.368e-2; 0.1946 -0.1257 -1.973e-4];
B = [0 0; -0.3295 0.304; -0.04073 -0.2737];
C = [-1.023 -1.444; 4.920 0.3648];
D = [-23.92; 5.921];

% Given
T = 120; % sec
u0 = 221; % ft/s

% Appendix, known
S = 5500; % ft^2
b = 198.68; % ft
c = 27.31; % ft
rho = 2.3769e-3; % lb sec^2/ft^4
W = 5.64e5; % lb
g = 32.2; % ft/s^2

%% Problem
omega = (2*pi)/T;
phi = atan2(omega*u0, g);
phiDeg = rad2deg(phi)
n = sec(phi);

CW = W/(0.5*rho*u0^2*S);

b1 = B*[0; -cos(phi)]*(omega*b/(2*u0));
x1 = A\b1;
beta = rad2deg(x1(1))
del_r = rad2deg(x1(2))
del_a = rad2deg(x1(3))

b2 = -D*(omega*c*sin(phi)/(2*u0)) + [0; (n-1)*CW];
x2 = C\b2;
delta_alpha = rad2deg(x2(1))
delta_del_e = rad2deg(x2(2))






