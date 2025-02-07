%% ASEN 3112 Lab 1 Problem 1 Script (Closed Thin Wall Analysis)
% Section 013 Group 6 - Liam Casey, Ian Faber, Chayce Johnson, Tori Morgheim, Ethan Stapp
% Code by: Tori Morgheim, Ian Faber

%% Housekeeping
clc; clear; close all;

%% Constants
dExt = 0.75; % Exterior diameter, inches
rExt = 3/8; % Exterior radius, inches
t = 1/16; % Wall thickness, inches
L = 1; % Extensometer gauge length, inches
G = 3.75 * 10^6; % Shear modulus, psi
griplength=10; %in

%% Data Extraction
CTWData = importdata('400inlbf_05.txt');

time = CTWData.data(:,1);
torsionalAngle = CTWData.data(:,2)*(pi/180);
shearStrain = CTWData.data(:,3);
torque = CTWData.data(:,4);
axialStrain = CTWData.data(:,5);
dphi=max(torsionalAngle)-min(torsionalAngle);
x_vec=axialStrain.*L+L;
dL=max(x_vec)-min(x_vec);
dphi_dx_ex=torsionalAngle./L;
dphi_dx_en=torsionalAngle./griplength;
%% Calculate shear strain from torsional Angle
% phi = gamma*L/Re
gamma = torsionalAngle*rExt/L;

%% Calculate Least Squares fit of both shear strain measurements


% [coef1, approxCurve1] = leastSquares(torque, shearStrain, 1);
% GJ_1 = coef1(1)*rExt;
% label1 = sprintf("Line of best fit: \n \\gamma = %.3e*T + %.3e", coef1(1), coef1(2));
% 
% [coef2, approxCurve2] = leastSquares(torque, gamma, 1);
% GJ_2 = coef2(1)*rExt;
% label2 = sprintf("Line of best fit: \n \\gamma = %.3e*T + %.3e", coef2(1), coef2(2));

[coef1, approxCurve1] = leastSquares(dphi_dx_ex, torque, 1);
GJ_1 = coef1(1);
[coef2, approxCurve2] = leastSquares(dphi_dx_en, torque, 1);
GJ_2 = coef2(1);


%% Plotting
figure
hold on;
title("OTW Extensometer Shear Strain vs. Torque")
plot(shearStrain, torque, 'b.', 'MarkerSize', 1);
%plot(shearStrain, approxCurve1(shearStrain), 'k--')
ylabel("Torque (in-lbf)")
xlabel("Shear Strain (rad)")
legend("Raw extensometer data",  'Location', 'best')

figure
hold on;
title("OTW Calculated Shear Strain vs. Torque")
plot((gamma) , torque, 'r.')
%plot((gamma), approxCurve2(gamma), 'k--')
ylabel("Torque (in-lbf)")
xlabel("Calculated Shear Strain (rad)")
legend("Calculated shear strain",  'Location', 'best')

%% Analytical calculations 
Ae=pi*rExt^2;
p=2*pi*rExt;
J_closed_analytical=(4*Ae^2*t)/p;
GJ_closed_an=G*J_closed_analytical;

Ri=rExt-t;
J_exact_analytical=(pi/2)*((rExt^4)-(Ri^4));
GJ_exact_an=G*J_exact_analytical;


%% Calculating Percent Error:
percentError1 = (abs((GJ_1-GJ_closed_an))/GJ_closed_an).*100;
percentError2 = (abs((GJ_2-GJ_exact_an))/GJ_exact_an).*100;

function [X,f] = leastSquares(t,y,p)
    % for writing this function, some skeleton code has been provided to
    % help you design the function to serve your purposes
    A = [];
    % write an expression for A, the input matrix
    for ii = 0:p
        col = t.^ii;
        A = [col, A];
    end
    % compute coefficient vector, x_hat
    x_hat = A\y;
    X = x_hat;
    
    % do not change the following lines of code. This will generate the
    % anonymous function handle "f" for you
    f = '@(x)';
    for i = 0:p
        f = strcat(f,'+',strcat(string(x_hat(i+1)),'.*x.^',string(p-i)));
    end
    eval(strcat('f = ',f,';'))
    
    while length(x_hat) < 7
        x_hat = [0;x_hat];
    end
    % workaround for MATLAB grader
%     f = @(x) x_hat(1)*x.^6 + x_hat(2)*x.^5 + x_hat(3)*x.^4 + x_hat(4)*x.^3 + x_hat(5)*x.^2 + x_hat(6)*x + x_hat(7);
    
end
