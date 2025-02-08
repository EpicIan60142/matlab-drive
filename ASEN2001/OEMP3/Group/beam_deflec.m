%% Beam deflect
% OEMP 3 Group 2 Part 5 a-c
% Group 11, 10:40 Lab
%
% Code by: Daniel Mascarenas
clear;
clc;


Shear_Max = 27*10^3; % Maximum Transverse sheat stress
BendMax = 35*10^3; % Maximum Bending stress

%Length in inches
B = 3.75; % Base  
H = 11.25; % Height
Area = 40.25 ; % (in^2)
FOS = 1.5153; % Ideal Factor of Saftey
L = 27.25 *12; % Wing length (in)
t = 1.75; % Thickness (in)

w0 = 2001; %Beam weight

p = 0.098; % lb /in^3

w1 = 0.5*L*(2001/12);% Maximum lift force
w2 = -p*Area; % Maximum weight force

E = 9.9*10^6;  % psi

I = ((B*H^3)/12)-(((H-2*t)*((B-2*t)^2))/12);% Moment of Inertia for hollow square


V = ((w0/2)*(-L))+(w2*(L^2)); % Takes place at x = 0
Q = (4.49)*(Area/2); % 4.49 = y'

Shearstress = (V*Q)/(I*t*2) % Min shear

M = ((w0/2)*((L^2)/3))+((w2/2)*(L^2)); %Min bending stress

sig = (-1*(M*(H/2))/I) %min bending stress

deflection = ((w0*L^4)/(30*E*I)) + ((w2*L^4)/(8*E*I)) %Beam Deflection in inches

Fos_bend = abs(sig/Shear_Max) % Min shear FOS
Fos_shear = abs(Shearstress/BendMax)%Min bending stress FOS

feet = deflection/12 %Deflection in feet

