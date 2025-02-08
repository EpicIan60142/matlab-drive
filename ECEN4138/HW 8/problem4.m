%% ECEN 4138 HW 8 Problem 5.27
% - Ian Faber

%% Housekeeping
clc; clear; close all

%% Setup vector of possibilities
numPoints = 100;
maxNum = 100;

p = linspace(0, maxNum, numPoints);
z = linspace(0, maxNum, numPoints);
% k = linspace(0, maxNum, numPoints);
k = 20;

%% Test conditions
epsilon = 1e-12;

[p, k, z] = meshgrid(p, k, z);

% p conditions
idxp1 = p > -(k/4) - 1 + sqrt(k.^2+8*k.*(z-1)+16)/4;
idxp2 = p < 0.1*k.*z;
idxp3 = p > -2;
idxp4 = (p <= z+epsilon) & (p >= z-epsilon);

% k conditions
idxk1 = k >= -4*(z-1) + 4*sqrt(z.*(z-2));
idxk2 = k > 0;
% idxk3 = k == 2;

% z conditions
idxz1 = z >= 2;
idxz2 = z > p;

% Compare all conditions
idxComp = idxp1.*idxp2.*idxp3.*idxk1.*idxk2.*idxz1.*idxz2;

idx = find(idxComp > 0);

pOpt = p(idx);
kOpt = k(idx);
zOpt = z(idx);

opt = [zOpt, pOpt, kOpt];

%% Make and test controller
s = tf('s');

G = 1/(s*(s+2));

vals = opt(7,:);

z = vals(1);
p = vals(2);
k = vals(3);
% k = 2;
% z = 2;
% p = 2;

C = k*(s+z)/(s+p);

L = minreal(C*G);

T = minreal(feedback(L,1))

roots(T.Denominator{:})

t = 0:0.001:10;

resp = lsim(T, t, t);
lsim(T,t,t);

resp(end) - t(end)

