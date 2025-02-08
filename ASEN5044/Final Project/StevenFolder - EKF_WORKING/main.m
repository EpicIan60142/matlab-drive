clc 
clear
close all

options = odeset('RelTol',1e-12, "AbsTol",1e-12);

%time vars
totalTime = 14000;
dT = 10;
thist = 0:10:14000;

%constants
mu = 398600;
rE = 6378;
omegaE = 2*pi/86400;
r0 = 6678;

%initial states
perturbation = [0;.075;0;-0.021];
x0 = [6678; 0;0; 6678*sqrt(mu/(6678^3))] + 0*perturbation;

%nominal solution
Xnom = r0*cos(sqrt(mu/r0^3).*thist);
Xdotnom = -sqrt(mu/r0^3)*r0*sin(sqrt(mu/r0^3).*thist);
Ynom = r0*sin(sqrt(mu/r0^3).*thist);
Ydotnom = r0*sqrt(mu/r0^3)*cos(sqrt(mu/r0^3).*thist);

nomSolution = [Xnom; Xdotnom; Ynom; Ydotnom];


% Assign input structures and matrices
load("orbitdeterm_finalproj_KFdata.mat")
P0 = diag([1,1e-3,1,1e-3]);

Gamma = [0 0; 1 0; 0 0; 0 1];


%constants
consts(1) = rE;
consts(2) = omegaE;
consts(3) = mu;

%LKF
% Q_LKF = 1e0*[
%                 5       1e-3    1e-3      0
%                 1e-3    1e-1    0        1e-3
%                 1e-3    0       5        1e-3
%                 0       1e-3    1e-3     1e-1
%              ];
Q_LKF = 1e-5*diag([1 1e-3 1 1e-3]);

filterParams = struct('R', Rtrue, 'Q', Q_LKF);
filterInit = struct('dx0', perturbation, 'P0', P0, 'xNom', nomSolution, "dt", dT);

LKF_init = LKF(ydata, filterParams, filterInit);

P0 = LKF_init.P_KF(:,:,5); % Initialize P0 with the 5th LKF estimate

%EKF
if true
    filterParams.Gamma = Gamma;
    filterParams.Q = 1e-5*diag([1 1e-3 1 1e-3]); %0.5*eye(4);
    filterParams.R = Rtrue;
    filterParams.dT = dT;
    filterParams.tvec = tvec;
    filterParams.consts = consts;
end

%initialize filter
perturb_x0 = [0;0.075;0;-0.021];
xhat0 = x0;


filterInit.xhat0 = xhat0;
filterInit.P0 = P0;

timesteps = 540;

% [epsilon_x_bar, epsilon_y_bar, NEESbounds, NISbounds] = chi2tests(@LKF, 50, 0.05, filterParams, Qtrue, Rtrue, filterInit, timesteps, x0);
[epsilon_x_bar, epsilon_y_bar, NEESbounds, NISbounds] = chi2tests(@EKF_Func, 100, 0.05, filterParams, Qtrue, Rtrue, filterInit, timesteps, x0);

% EKFdata = EKF_Func(ydata, filterParams, filterInit);

% Check alpha
NEES_consistency = epsilon_x_bar >= NEESbounds(1) & epsilon_x_bar <= NEESbounds(2);
NEES_alpha = 1 - sum(NEES_consistency)/timesteps
NIS_consistency_1_stat = epsilon_y_bar(1,:) >= NISbounds(1,1) & epsilon_y_bar(1,:) <= NISbounds(1,2);
NIS_alpha_1_stat = 1 - sum(NIS_consistency_1_stat)/timesteps
NIS_consistency_2_stat = epsilon_y_bar(2,:) >= NISbounds(2,1) & epsilon_y_bar(2,:) <= NISbounds(2,2);
NIS_alpha_2_stat = 1 - sum(NIS_consistency_2_stat)/timesteps

t = 1:(timesteps + 1);

figure()
hold on
title("NIS - 1 station")
plot(t,epsilon_y_bar(1,:), "bo")
yline(NISbounds(1,1), "k--")
yline(NISbounds(2,1), "k--")
ylim([0,10])

figure()
hold on
title("NIS - 2 stations")
plot(epsilon_y_bar(2,:), "bo")
yline(NISbounds(1,2), "k--")
yline(NISbounds(2,2), "k--")
ylim([0,10])

figure()
hold on
title("NEES")
plot(1:timesteps+1,epsilon_x_bar, "ro")
yline(NEESbounds(1), "k--")
yline(NEESbounds(2), "k--")







