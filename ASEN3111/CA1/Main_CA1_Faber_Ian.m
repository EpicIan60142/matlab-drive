%% ASEN 3111 - Computational Assignment 01 - Main
%   Driver program for CA 01 - Computation of Lift and Drag
%       Program that calculates coefficients of lift and drag for 2 shapes:
%       A rotating cylinder and a NACA 0012 airfoil.
%
%   Author: Ian Faber
%   Collaborators: Maggie Wussow, Katy McCutchan, Anna Sophia Rorrer
%   Warren, Jonathan Abrams
%   Date Started: 1/21/2023
%   Date Finished: 2/4/2023
%

%% Housekeeping
clc; clear; close all;

%% Starter Command Window Management
fprintf("\t+--------------------------------+\n")
fprintf("\t\tCA 1 Results - Ian Faber\n\n\n")
fprintf("\tProblem 1: \n\n")

%% Problem 1
% Rotating Cylinder drag and lift

% Constants
nRef = 1000;
nPanels = 10;
clTrapz = [];
cdTrapz = [];
clSimp = [];
cdSimp = [];

% Reference Cp curve
thetaRef = linspace(0, 2*pi, nRef+1);
cpRef = -4*sin(thetaRef).*(sin(thetaRef)+1);
clRef = cpRef.*sin(thetaRef);
cdRef = cpRef.*cos(thetaRef);

cp = @(theta) -4*sin(theta).*(sin(theta)+1);
cl = @(theta) cp(theta).*sin(theta);
cd = @(theta) cp(theta).*cos(theta);

% Analytical coefficient of lift
    % Calculated by hand, c_l = -0.5*int{0}^{2*pi}(c_p(Theta)*sin(Theta))dTheta
clAnalytical = 2*pi;
fprintf("Analytical coefficient of lift: 2*pi or %.6f \n", clAnalytical)

% Analytical coefficient of drag
    % Calculated by hand, c_d = -0.5*int{0}^{2*pi}(c_p(Theta)*cos(Theta))dTheta
cdAnalytical = 0;
fprintf("Analytical coefficient of drag: 0 or %.3f \n", cdAnalytical)

% Composite Trapezoidal
for k = 1:nPanels
    theta = linspace(0, 2*pi, k + 1);
    clTrapz = [clTrapz; -0.5*compTrapz(theta, cl)];
    cdTrapz = [cdTrapz; -0.5*compTrapz(theta, cd)];
end

clTrapzErr = (clTrapz - clAnalytical)/clAnalytical;

% Composite Simpson's
for k = 1:nPanels
    theta = linspace(0, 2*pi, 2*k + 1);
    clSimp = [clSimp; -0.5*compSimp(theta, cl)];
    cdSimp = [cdSimp; -0.5*compSimp(theta, cd)];
end

clSimpErr = (clSimp - clAnalytical)/clAnalytical;

% Number of panels for 1% trapezoidal error
for k = 1:length(clTrapzErr)
    if abs(clTrapzErr(k)) < 0.01
        nTrapz = k;
        break;
    end
end

% nTrapz
fprintf("Number of panels needed for 1%% relative trapezoidal error: %.0f \n", nTrapz)

% Number of panels for 1% simpson's error
for k = 1:length(clSimpErr)
    if abs(clSimpErr(k)) < 0.01
        nSimp = k;
        break;
    end
end

% nSimp
fprintf("Number of panels needed for 1%% relative simpson's error: %.0f \n", nSimp)

% Plotting

% figure
% hold on
% title("C_p")
% plot(thetaRef, cpRef)
% plot(theta, cp(theta), 'k--')
% xlabel("Theta (rad)");
% ylabel("C_p");
% legend("Reference", "Trapezoidal", 'Location', 'best')
% 
% figure
% hold on
% title("$C_p\sin(\theta)$", "Interpreter","latex")
% plot(thetaRef, clRef)
% plot(theta, cl(theta))
% xlabel("Theta (rad)");
% ylabel("C_l");
% legend("Reference", "Trapezoidal", 'Location', 'best')
% 
% figure
% hold on
% title("$C_p\cos(\theta)$", "Interpreter","latex")
% plot(thetaRef, cdRef)
% plot(theta, cd(theta))
% xlabel("Theta (rad)");
% ylabel("C_d");
% legend("Reference", "Trapezoidal", 'Location', 'best')

figure
hold on
title("Composite Trapezoidal c_l");
plot(clTrapz);
yline(clAnalytical, 'k--');
xlabel("Number of panels");
ylabel("c_l");
legend("Trapezoidal", "Analytical", 'Location', 'best')

figure
hold on
title("Composite Trapezoidal c_d");
plot(cdTrapz);
yline(cdAnalytical, 'k--');
xlabel("Number of panels")
ylabel("c_d")
legend("Trapezoidal", "Analytical", 'Location', 'best')

figure
hold on
title("Composite Simpson's c_l");
plot(clSimp);
yline(clAnalytical, 'k--');
xlabel("Number of panels");
ylabel("c_l");
legend("Simpson's", "Analytical", 'Location', 'best')

figure
hold on
title("Composite Simpson's c_d");
plot(cdSimp);
yline(cdAnalytical, 'k--');
xlabel("Number of panels")
ylabel("c_d")
legend("Simpson's", "Analytical", 'Location', 'best')

%% Intermediate Command Window management
fprintf("\n\n\tProblem 2: \n\n")

%% Problem 2

% Turbulent flow over NACA 0012 airfoil

% Constants
alpha = deg2rad(10); % deg
c = 3; % m
V_inf = 50; % m/s
rho_inf = 1.225; % kg/m^3
p_inf = 101.3e3; % Pa
Re = 10e6;
xRef = linspace(0, c, 1000);
t = 0.12; % xx/100

load Cp

nPointsMax = 5000;
nTrue = 1000000;

cn_up = zeros(nPointsMax, 1);
cn_low = zeros(nPointsMax, 1);
ca_up = zeros(nPointsMax, 1);
ca_low = zeros(nPointsMax, 1);
cl = zeros(nPointsMax, 1);
cd = zeros(nPointsMax, 1);

% Airfoil curve
y = @(x, t) (t/0.2)*c*(0.2969*sqrt(x/c) - 0.126*(x/c) - 0.3516*(x/c).^2 + 0.2843*(x/c).^3 -0.1036*(x/c).^4);
yRef = y(xRef, t);


for k = 1:nPointsMax
    xc = linspace(0, 1, k);
    x = linspace(0,c,k);
    cp_up = fnval(Cp_upper, xc);
    cp_low = fnval(Cp_lower, xc);
    cn_up(k) = compTrapz(x, @(x) -cp_up);
    cn_low(k) = compTrapz(x, @(x) cp_low);
    cn = (1/c)*(cn_up + cn_low);
    
    ca_up(k) = compTrapz(y(x,t), @(x) cp_up);
    ca_low(k) = compTrapz(-y(x,t), @(x) -cp_low);
    ca = (1/c)*(ca_up + ca_low);

    cl = cn*cos(alpha) - ca*sin(alpha);
    cd = cn*sin(alpha) + ca*cos(alpha);
end

L = 0.5*cl*rho_inf*c*(V_inf)^2;
D = 0.5*cd*rho_inf*c*(V_inf)^2;

% "True" values
xTrue = linspace(0,c,nTrue);
xcTrue = linspace(0, 1, nTrue);
cp_upTrue = fnval(Cp_upper, xcTrue);
cp_lowTrue = fnval(Cp_lower, xcTrue);
cn_upTrue = compTrapz(xTrue, @(xTrue) -cp_upTrue);
cn_lowTrue = compTrapz(xTrue, @(xTrue) cp_lowTrue);
cnTrue = (1/c)*(cn_upTrue + cn_lowTrue);

ca_upTrue = compTrapz(y(xTrue,t), @(xTrue) cp_upTrue);
ca_lowTrue = compTrapz(-y(xTrue,t), @(xTrue) -cp_lowTrue);
caTrue = (1/c)*(ca_upTrue + ca_lowTrue);

clTrue = cnTrue*cos(alpha) - caTrue*sin(alpha);
cdTrue = cnTrue*sin(alpha) + caTrue*cos(alpha);

LTrue = 0.5*clTrue*c*rho_inf*(V_inf)^2;
DTrue = 0.5*cdTrue*c*rho_inf*(V_inf)^2;

for k = 1:length(L)
    errLift = abs((L(k) - LTrue)/LTrue);
    if errLift < 0.01
        nLift = k;
        break;
    end
end

first = false;
for k = 1:length(D)
    errDrag = abs((D(k) - DTrue)/DTrue);
    if errDrag < 0.01
        if first
            nDrag = k;
            break;
        end
        first = true;
    end
end

fprintf("Lift solution: %.3f N\n", LTrue);
fprintf("Number of integration points needed for 1%% lift error: %.0f \n", 2*nLift - 2); % Avoid double counting LE and TE
fprintf("Drag solution: %.3f N\n", DTrue);
fprintf("Number of integration points needed for 1%% drag error: %.0f \n", 2*nDrag - 2); % Avoid double counting LE and TE

% Plotting

% % Upper curve
% for k = 1:nPointsMax
%     xc = linspace(0, 1, k);
%     cp_up = fnval(Cp_upper, xc);
%     cp = @(x) cp_up;
%     clTrapzUp = [clTrapzUp; compTrapz(xc, cp)];
% end
% 
% % Lower curve
% for k = 1:nPointsMax
%     xc = linspace(0, 1, k);
%     cp_low = fnval(Cp_lower, xc);
%     cp = @(x) cp_low;
%     clTrapzLow = [clTrapzLow; compTrapz(xc, cp)];
% end
% 
% clTrapz2 = clTrapzLow - clTrapzUp;

x = linspace(0,2*nPointsMax - 2,length(L));

figure
hold on
title("Lift Solution vs. Number of Integration Points")
plot(x,L)
xlabel("Number of integration points")
ylabel("Lift (N)")

figure
hold on
title("Drag Solution vs. Number of Integration Points")
plot(x,D)
xlabel("Number of integration points")
ylabel("Drag (N)")

% figure
% hold on;
% title("NACA 0012 Airfoil")
% plot(x, -cp_up, 'r')
% plot(x, -cp_low, 'm')
% plot(xRef, yRef, 'b-');
% plot(xRef, -yRef, 'b-');
% ylim([-c/2, 2*c])
% xlabel("x")
% legend("Upper Cp", "Lower Cp", "Airfoil")

%% Final Command Window Management
fprintf("\n\n\t+--------------------------------+\n\n")
disp(getRalphie)