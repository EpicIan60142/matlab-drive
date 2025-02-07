%% ASEN 3111 - Computational Assignment 03 - Main
%   Driver program for CA 03 - Flow over Thick Airfoils and Finite Wings
%       Program that explores the effectiveness of thin airfoil theory and
%       compares it to the vortex panel method. Also explores PLLT (Prandtl
%       Lifing Line Theory) and its applications.
%
%   Author: Ian Faber
%   Collaborators: Maggie Wussow
%   Date started: 03/12/2023
%   Data finished: 03/14/2023

%% Housekeeping
clc; clear; close all;

tic;

%% Starter Command Window Management
fprintf("\t+--------------------------------+\n")
fprintf("\t\t CA 3 Results - Ian Faber\n\n")

%% Problem 1
fprintf("\n\t Problem 1:\n\n")

V_inf = 200; % ft/s
alpha = 10; % deg
c = 3.2; % ft

% lower = 2;
% upper = 1;
%
% [x,y,xVort,yVort,~] = NACA_Airfoils(0,0,12,c,5000);
% 
% clTrue = Vortex_Panel(xVort, yVort, V_inf, alpha)
%
% save("Data_CA3_Faber_Ian.mat", "clTrue");
%
% xLower = x(lower,:);
% xUpper = x(upper,:);
% 
% yLower = y(lower,:);
% yUpper = y(upper,:);
% 
% figure
% hold on
% axis equal
% plot(xLower,yLower)
% plot(xUpper,yUpper)

load("Data_CA3_Faber_Ian.mat")

cl = 999999;
N = 4;
while abs((cl-clTrue)/clTrue) > 0.01
    N = N + 1;

    [x,y,xVort,yVort,~] = NACA_Airfoils(0,0,12,c,N);
    cl = Vortex_Panel(xVort, yVort, V_inf, alpha);

end

fprintf("Number of panels needed for 1%% relative error: %.0f\n", N);

%% Problem 2
fprintf("\n\tProblem 2:\n\n")

alpha = -10:1:10;
cl0012_TAT = 2*pi*deg2rad(alpha);

[~,~,xVort0006,yVort0006,~] = NACA_Airfoils(0,0,06,c,N);
[~,~,xVort0012,yVort0012,~] = NACA_Airfoils(0,0,12,c,N);
[~,~,xVort0024,yVort0024,~] = NACA_Airfoils(0,0,24,c,N);

for k = 1:length(alpha)
    cl0006(k) = Vortex_Panel(xVort0006, yVort0006, V_inf, alpha(k));
    cl0012(k) = Vortex_Panel(xVort0012, yVort0012, V_inf, alpha(k));
    cl0024(k) = Vortex_Panel(xVort0024, yVort0024, V_inf, alpha(k));
end

[coef0012_TAT,~] = polyfit(alpha,cl0012_TAT,1);
[coef0006,~] = polyfit(alpha,cl0006,1);
[coef0012,~] = polyfit(alpha,cl0012,1);
[coef0024,~] = polyfit(alpha,cl0024,1);

fprintf("TAT lift-curve slope: %.4f/deg, TAT zero lift AoA: %.4f\n", coef0012_TAT)
fprintf("NACA 0006 lift-curve slope: %.4f/deg, NACA 0006 zero lift AoA: %.4f deg\n", coef0006)
fprintf("NACA 0012 lift-curve slope: %.4f/deg, NACA 0012 zero lift AoA: %.4f deg\n", coef0012)
fprintf("NACA 0024 lift-curve slope: %.4f/deg, NACA 0024 zero lift AoA: %.4f deg\n", coef0024)

figure
hold on
grid on
title("Lift-Curve Slope for Symmetric Airfoils")
plot(alpha, cl0012_TAT, 'k--')
plot(alpha, cl0006);
plot(alpha, cl0012);
plot(alpha, cl0024);
xline(0,'k-')
yline(0,'k-')
xlabel("\alpha (deg)")
ylabel("c_l")

legend("TAT","NACA 0006", "NACA 0012", "NACA 0024", 'Location', 'best')


%% Problem 3
fprintf("\n\tProblem 3:\n\n")

% Vortex Panel
[~,~,xVort2412,yVort2412,dzdx2412] = NACA_Airfoils(2,4,12,c,N);
[~,~,xVort4412,yVort4412,dzdx4412] = NACA_Airfoils(4,4,12,c,N);

for k = 1:length(alpha)
    cl2412(k) = Vortex_Panel(xVort2412, yVort2412, V_inf, alpha(k));
    cl4412(k) = Vortex_Panel(xVort4412, yVort4412, V_inf, alpha(k));
end

[coef2412,~] = polyfit(alpha,cl2412,1);
[coef4412,~] = polyfit(alpha,cl4412,1);

% TAT
theta = linspace(0,pi,N/2 + 1);
cl2412_TAT = 2*pi*(deg2rad(alpha) - (1/pi)*trapz(theta,dzdx2412(theta).*(1-cos(theta))));
cl4412_TAT = 2*pi*(deg2rad(alpha) - (1/pi)*trapz(theta,dzdx4412(theta).*(1-cos(theta))));

[coef2412_TAT,~] = polyfit(alpha,cl2412_TAT,1);
[coef4412_TAT,~] = polyfit(alpha,cl4412_TAT,1);

alphaL0 = @(coef) -coef(2)/coef(1); % Solve for x-intercept

fprintf("NACA 0012:\n")
fprintf("\tTAT: lift-curve slope: %.4f/deg, zero-lift AoA: %.4f deg\n", coef0012_TAT(1), alphaL0(coef0012_TAT));
fprintf("\tVPM: lift-curve slope: %.4f/deg, zero-lift AoA: %.4f deg\n\n", coef0012(1), alphaL0(coef0012));
fprintf("NACA 2412:\n")
fprintf("\tTAT: lift-curve slope: %.4f/deg, zero-lift AoA: %.4f deg\n", coef2412_TAT(1), alphaL0(coef2412_TAT));
fprintf("\tVPM: lift-curve slope: %.4f/deg, zero-lift AoA: %.4f deg\n\n", coef2412(1), alphaL0(coef2412));
fprintf("NACA 4412:\n")
fprintf("\tTAT: lift-curve slope: %.4f/deg, zero-lift AoA: %.4f deg\n", coef4412_TAT(1), alphaL0(coef4412_TAT));
fprintf("\tVPM: lift-curve slope: %.4f/deg, zero-lift AoA: %.4f deg\n", coef4412(1), alphaL0(coef4412));

figure
hold on
grid on
title("Lift-Curve Slope for Cambered Airfoils")
a = plot(alpha, cl0012, 'b-');
plot(alpha, cl0012_TAT ,'b--')
b = plot(alpha, cl2412, 'r-');
plot(alpha, cl2412_TAT, 'r--')
c = plot(alpha, cl4412, 'k-');
d = plot(alpha, cl4412_TAT, 'k--');
xline(0,'k-')
yline(0,'k-')
xlabel("\alpha (deg)")
ylabel("c_l")

subset = [a,b,c,d];

legend(subset, ["NACA 0012", "NACA 2412", "NACA 4412", "TAT; All dashed lines"], 'Location', 'best');

%% Problem 4
fprintf("\n\tProblem 4:\n\n")

alpha = 0;

geo_r = 0 + alpha; % deg
geo_t = 0 + alpha; % deg
aero_r = rad2deg(alphaL0(coef0012)); % deg
aero_t = rad2deg(alphaL0(coef0012)); % deg
a0_r = rad2deg(coef0012(1)); % 1/rad, had to invert the deg2rad function since 1/deg = (pi/180) rad / deg
a0_t = rad2deg(coef0012(1)); % 1/rad

AR = [4,6,8,10]; % Aspect ratios of interest
ctcr = 0:0.01:1; % Vector of ratios to test
c_r = 1; % ft, fixed root chord
N = 50; % Number of odd terms

figure
hold on
title("$\delta$ vs. $\frac{c_t}{c_r}$",'Interpreter','latex')

for kk = 1:length(AR)
    for k = 1:length(ctcr)
        c_t = ctcr(k)*c_r; % ft
        c = @(theta) c_r + (c_t - c_r)*cos(theta); % ft
        
        b = 0.5*AR(kk)*(c_r + c_t);
    
        [e, ~, ~] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N);
        delta(k) = (1/e)-1;
    end

    plot(ctcr,delta)
    label(kk) = sprintf("AR = %.0f", AR(kk));
end

xlabel("$\frac{c_t}{c_r}$",'Interpreter','latex')
ylabel("\delta")
legend(label, 'Location', 'best')

fprintf("See plot for reproduction of Figure 5.20 from Anderson\n")


%% Problem 5
fprintf("\n\tProblem 5:\n\n")

alpha = 3; % deg
V_inf = 60*1.68781; % knots -> ft/s
h = 15000; % ft

[~,~,~,rho] = atmosisa(h*0.3048); % ft -> m
rho = rho*0.00194032; % kg/m^3 -> slug/ft^3

b = 33.333; % ft
c_r = 5.333; % ft
c_t = 3.708333; % ft
geo_r = 1 + alpha; % deg
geo_t = 0 + alpha; % deg
aero_r = alphaL0(coef2412); % deg
aero_t = alphaL0(coef0012); % deg
a0_r = coef2412(1)*(180/pi); % 1/rad, had to invert the deg2rad function since 1/rad = 1 / (deg * pi/180 rad/deg) = 180/pi * 1/deg
a0_t = coef0012(1)*(180/pi); % 1/rad, ditto

c = @(theta) c_r + (c_t - c_r)*cos(theta);
wing = linspace(0,pi/2,1000);
S = 2*trapz(wing, c(wing).*(b/2).*sin(wing));

q_inf = 0.5*rho*V_inf^2;

% NTrue = 10000;
% [~, c_LTrue, c_DiTrue] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, NTrue);
% save("Data_CA3_Faber_Ian.mat","c_LTrue","c_DiTrue",'-append')

LTrue = q_inf*c_LTrue*S;
DiTrue = q_inf*c_DiTrue*S;

fprintf("Lift and Induced Drag solution:\n")
fprintf("\tLift solution: %.3f lb\n", LTrue);
fprintf("\tInduced Drag solution: %.3f lb\n\n", DiTrue);

c_L = 99999;
c_Di = 99999;
N = 0;

% Detection booleans
tenPercentCL = false;
onePercentCL = false;
pOnePercentCL = false;
tenPercentCD = false;
onePercentCD = false;
pOnePercentCD = false;

while (abs((c_L - c_LTrue)/c_LTrue) > 0.001) || (abs((c_Di - c_DiTrue)/c_DiTrue) > 0.001)
    N = N + 1;
    
    [~, c_L, c_Di] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N);

    % Check C_L convergence
    if abs((c_L - c_LTrue)/c_LTrue) < 0.1 && ~tenPercentCL
        NTenPercentCL = N;
        tenPercentCL = true;
    end

    if abs((c_L - c_LTrue)/c_LTrue) < 0.01 && ~onePercentCL
        NOnePercentCL = N;
        onePercentCL = true;
    end

    if abs((c_L - c_LTrue)/c_LTrue) < 0.001 && ~pOnePercentCL
        NpOnePercentCL = N;
        pOnePercentCL = true;
    end

    % Check C_Di convergence
    if abs((c_Di - c_DiTrue)/c_DiTrue) < 0.1 && ~tenPercentCD
        NTenPercentCD = N;
        tenPercentCD = true;
    end

    if abs((c_Di - c_DiTrue)/c_DiTrue) < 0.01 && ~onePercentCD
        NOnePercentCD = N;
        onePercentCD = true;
    end

    if abs((c_Di - c_DiTrue)/c_DiTrue) < 0.001 && ~pOnePercentCD
        NpOnePercentCD = N;
        pOnePercentCD = true;
    end
end

fprintf("Lift:\n")
fprintf("\tNumber of odd terms needed for 10%% Lift error: %.0f\n", NTenPercentCL);
fprintf("\tNumber of odd terms needed for 1%% Lift error: %.0f\n", NOnePercentCL);
fprintf("\tNumber of odd terms needed for 0.1%% Lift error: %.0f\n\n", NpOnePercentCL);

fprintf("Drag:\n")
fprintf("\tNumber of odd terms needed for 10%% Induced Drag error: %.0f\n", NTenPercentCD);
fprintf("\tNumber of odd terms needed for 1%% Induced Drag error: %.0f\n", NOnePercentCD);
fprintf("\tNumber of odd terms needed for 0.1%% Induced Drag error: %.0f\n\n", NpOnePercentCD);

%% End Code
toc;

%% Final Command Window Management
fprintf("\n\n\t+--------------------------------+\n\n")
disp(getRalphie)