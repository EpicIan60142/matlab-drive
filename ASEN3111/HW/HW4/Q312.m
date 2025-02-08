%% ASEN 3111 HW 4 Q3.12
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Code
V_inf = 5;
lambda = 2*pi*V_inf; % stagnation 1 ft upstream
thres = 50;

nPoints = 100;
xVec = linspace(-2, 2, nPoints);
yVec = linspace(-3, 3, nPoints);

[x,y] = meshgrid(xVec, yVec);

Psi = V_inf*y + (lambda/(2*pi))*mod(atan2(y,x),2*pi);
Phi = V_inf*x + (lambda/(2*pi))*log(sqrt(x.^2 + y.^2));

[u, v] = gradient(Phi, mean(diff(xVec)), mean(diff(yVec)));
v(v>thres) = thres;
v(v<-thres) = -thres;
u(u>thres) = thres;
u(u<-thres) = -thres;
V = sqrt(u.^2 + v.^2);
Cp = 1 - (V/V_inf).^2;

levmin = min(Psi,[],'all');
levmax = max(Psi,[],'all');
levels = linspace(levmin,levmax,50)';

figure
hold on
title("Cp Distribution")
contour(x,y,Psi,levels,'LineWidth',1.5)
contour(x,y,Psi,1,'LineWidth',3,'LineColor','k')
colorbar

levmin = min(min(Cp));
levmax = max(max(Cp));
levels = linspace(levmin,levmax,50)';

figure
hold on
title("Cp Distribution")
contour(x,y,Cp,levels,'LineWidth',1.5)
contour(x,y,Cp,1,'LineWidth',3,'LineColor','k')
colorbar
plot(xVec, -Cp(25,:))
xlabel("x location")
ylabel("-Cp")