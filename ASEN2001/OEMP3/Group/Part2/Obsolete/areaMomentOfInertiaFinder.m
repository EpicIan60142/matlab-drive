%% OEMP 3 Group 2 Parts 2-3
%   Moment of inertia finder

clc; clear; close all;

%% Constants

w0 = 2001;
L = 27.25*12; % Convert to inches
rho = 0.16; % Titanium, lb/in^3
sigmaYield = 120000; % Titanium, 120 ksi

%% Equations

moment = @(A, x) (w0/2)*((-x^3)/(3*L) + x^2 - L*x + (L^2)/3) - (rho*A/2)*(L-x)^2;

maxPoint = @(A) (-L*(2*A - w0))/w0;

inertiaCircle = @(r) (pi*r^4)/4;
inertiaRectangle = @(b, h) (b*h^3)/12;

bendingStress = @(M, y, I) M*y/I;

factorOfSafety = @(sigmaYield, sigmaApplied) sigmaYield/sigmaApplied;

%% Circle

factorsCircle = [];
factorOneCircle = zeros(1,48);
for k = 1:48 % Go from 1/4 in to 12 in
    
    r = k/4;

    inertia = inertiaCircle(r);
    area = pi*r^2;
    maxX = maxPoint(area);
    sigmaBend = bendingStress(moment(area, maxX), r, inertia);

    factorOneCircle(k) = factorOfSafety(sigmaYield, sigmaBend);
end
factorsCircle = [factorsCircle, factorOneCircle];

factorsCircle







