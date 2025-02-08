% AUTOGRADER WILL CHECK FOR VARIABLES OF THESE NAMES
% Section 011 - Ian Faber, Kevin McGough, Alex Putman, Blake Wilson

% ONLY LOOK AT DATA FOR PROBLEM 3
% Derive IC's for each problem

m = 0.068; % kg
g = 9.81; % m/s^2
d = 0.06; % m
Km = 0.0024; % N*m/N
Ix = 6.85*10^-5; % kgm^2
Iy = 9.2*10^-5; % kgm^2
Iz = 1.35*10^-4; % kgm^2
nu = 10^-3; % N/(m/s)^2
mu = 2*10^-6; % M*m/(rad/s)^2 

% Problem 1
p0 = [0; 0; -5];
a0 = [0; 0; 0];
v0 = [0; 0; 0];
r0 = [0; 0; 0];
x0_steady_hover = [p0; a0; v0; r0];  %12x1 column vector of initial state for Problems 1 and 2a
Zc_hover = -m*g;%scalar body-z thrust for steady hover
forces_steady_hover = [Zc_hover/4; Zc_hover/4; Zc_hover/4; Zc_hover/4] ;%4x1 column vector of motor forces for steady hover

% Problem 2b
phi_b = deg2rad(2.1462559951889); %scalar euler angle phi for Problem 2b
Zc_b = -0.667548283334; %scalar body-z thrust for Problem 2b
v_E = 4.9964924247073;
w_E = -0.1827252369457;
p0 = [0; 0; -5];
a0 = [phi_b; 0; 0];
v0 = [0; v_E; w_E];
r0 = [0; 0; 0];
x0_2b = [p0; a0; v0; r0]; %12x1 column vector of initial state for Problem 2b
forces_2b = [Zc_b/4; Zc_b/4; Zc_b/4; Zc_b/4]; %4x1 column vector of motor forces for 2b

% Problem 2c
theta_c = -phi_b; %scalar euler angle theta for Problem 2c
Zc_c = Zc_b; %scalar body-z thrust for Problem 2c
u_E = v_E;
w_E_2 = w_E;
p0 = [0; 0; -5];
a0 = [0; theta_c; 0];
v0 = [u_E; 0; w_E_2];
r0 = [0; 0; 0];
x0_2c = [p0; a0; v0; r0]; %12x1 column vector of initial state for Problem 2c
forces_2c = [Zc_c/4; Zc_c/4; Zc_c/4; Zc_c/4]; %4x1 column vector of motor forces for 2c

