%% Housekeeping
clc 
clear
close all

%% Setup
options = odeset('RelTol',1e-12, "AbsTol",1e-12);

%time vars
totalTime = 14000;
dT = 10;
thist = 0:dT:14000;

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
P0 = diag([0.5,1e-3,0.5,1e-3]);

Gamma = [0 0; 1 0; 0 0; 0 1];


%constants
consts(1) = rE;
consts(2) = omegaE;
consts(3) = mu;

%LKF
Q_LKF_init = 5e-3*diag([1 1e-3 1 1e-3]);
Q_LKF = 5e-4*diag([1 1e-5 1 1e-5]);
Q_LKF_old = 1e-5*diag([1 1e-3 1 1e-3]);
Q_EKF = 1e-5*diag([1 1e-3 1 1e-3]);

filterParams = struct('R', Rtrue, 'Q', Q_LKF_init);
filterInit = struct('dx0', 0*perturbation, 'P0', P0, 'xNom', nomSolution, "dt", dT);

LKF_init = LKF(ydata, filterParams, filterInit);

P0_old = LKF_init.P_KF(:,:,5);
P0 = 1e-3*LKF_init.P_KF(:,:,5); % Initialize P0 with a scaled 5th LKF estimate

if false
    filterParams.Gamma = Gamma;
    filterParams.Q = Q_EKF; %0.5*eye(4);
    filterParams.R = Rtrue;
    filterParams.dT = dT;
    filterParams.tvec = tvec;
    filterParams.consts = consts;
else
    filterParams.Q = Q_LKF;
end

%% Run LKF
%initialize filter
perturb_x0 = [0;0.075;0;-0.021];
xhat0 = x0;

filterInit.xhat0 = xhat0;
filterInit.P0 = P0;

timesteps = 540;

[epsilon_x_bar, epsilon_y_bar, NEESbounds, NISbounds, xMonteCarlo, yMonteCarlo, filterData] = chi2tests(@LKF, 100, 0.05, filterParams, Qtrue, Rtrue, filterInit, timesteps, x0);
% [epsilon_x_bar, epsilon_y_bar, NEESbounds, NISbounds] = chi2tests(@EKF_Func, 100, 0.05, filterParams, Qtrue, Rtrue, filterInit, timesteps, x0);

% EKFdata = EKF_Func(ydata, filterParams, filterInit);
% LKFdata = LKF(ydata, filterParams, filterInit);

% filterData = LKFdata;

% Check alpha
NEES_consistency = epsilon_x_bar >= NEESbounds(1) & epsilon_x_bar <= NEESbounds(2);
NEES_alpha = 1 - sum(NEES_consistency)/timesteps
NIS_consistency_1_stat = epsilon_y_bar(1,:) >= NISbounds(1,1) & epsilon_y_bar(1,:) <= NISbounds(1,2);
NIS_alpha_1_stat = 1 - sum(NIS_consistency_1_stat)/timesteps
NIS_consistency_2_stat = epsilon_y_bar(2,:) >= NISbounds(2,1) & epsilon_y_bar(2,:) <= NISbounds(2,2);
NIS_alpha_2_stat = 1 - sum(NIS_consistency_2_stat)/timesteps

%% Plot NEES/NIS
t = 1:(timesteps + 1);

figure()
hold on
title("NIS - 1 station")
plot(t,epsilon_y_bar(1,:), "b.")
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

%% Plot noisy states
figure
sgtitle("Typical Monte Carlo Simulation States - Last Monte Carlo run")
ax(1) = subplot(2,2,1);
    hold on; grid on;
    title("X state")
    plot(t, xMonteCarlo(1,:))
    xlabel("Time [sec]"); ylabel("X [km]")
ax(2) = subplot(2,2,2);
    hold on; grid on;
    title("Xdot state")
    plot(t, xMonteCarlo(2,:))
    xlabel("Time [sec]"); ylabel("Xdot [km/s]")
ax(3) = subplot(2,2,3);
    hold on; grid on;
    title("Y state")
    plot(t, xMonteCarlo(3,:))
    xlabel("Time [sec]"); ylabel("Y [km]")
ax(4) = subplot(2,2,4);
    hold on; grid on;
    title("Ydot state")
    plot(t, xMonteCarlo(4,:))
    xlabel("Time [sec]"); ylabel("Ydot [km/s]")
linkaxes(ax, 'x');

%% Plot noisy data
kTime = 1:length(filterData.yMeas_KF);

figure
sgtitle("Typical Monte Carlo Simulation Measurements - Last Monte Carlo run")
ax(1) = subplot(3,1,1);
hold on; grid on;
title("\rho vs. time")
for k = 1:length(filterData.yMeas_KF)
    if ~isempty(filterData.yMeas_KF{k})
        plot(kTime(k), filterData.yMeas_KF{k}(1), 'b.')
        if length(filterData.yMeas_KF{k}) > 3
            plot(kTime(k), filterData.yMeas_KF{k}(4), 'b.')
        end
    end
end

ax(2) = subplot(3,1,2);
hold on; grid on;
title("\rho_{Dot} vs. time")
for k = 1:length(filterData.yMeas_KF)
    if ~isempty(filterData.yMeas_KF{k})
        plot(kTime(k), filterData.yMeas_KF{k}(2), 'b.')
        if length(filterData.yMeas_KF{k}) > 3
            plot(kTime(k), filterData.yMeas_KF{k}(5), 'b.')
        end
    end
end

ax(3) = subplot(3,1,3);
hold on; grid on;
title("\phi vs. time")
for k = 1:length(filterData.yMeas_KF)
    if ~isempty(filterData.yMeas_KF{k})
        plot(kTime(k), filterData.yMeas_KF{k}(3), 'b.')
        if length(filterData.yMeas_KF{k}) > 3
            plot(kTime(k), filterData.yMeas_KF{k}(6), 'b.')
        end
    end
end

linkaxes(ax, 'x')

%% Plot state errors
t_KF = filterData.t_KF;
x_KF = filterData.x_KF;
sig_KF = filterData.sig_KF;

stateError = nomSolution(:,1:size(x_KF,2)) - x_KF;

figure
sgtitle("LKF State Errors - Last Monte Carlo run")
ax(1) = subplot(2,2,1);
    hold on; grid on;
    title("X state error")
    plot(t_KF, stateError(1,:))
    plot(t_KF, 2*sig_KF(1,:), 'r--');
    plot(t_KF, -2*sig_KF(1,:), 'r--')
    xlabel("Time [sec]"); ylabel("X [km]")
ax(2) = subplot(2,2,2);
    hold on; grid on;
    title("Xdot state error")
    plot(t_KF, stateError(2,:))
    plot(t_KF, 2*sig_KF(2,:), 'r--');
    plot(t_KF, -2*sig_KF(2,:), 'r--')
    xlabel("Time [sec]"); ylabel("Xdot [km/s]")
ax(3) = subplot(2,2,3);
    hold on; grid on;
    title("Y state error")
    plot(t_KF, stateError(3,:))
    plot(t_KF, 2*sig_KF(3,:), 'r--');
    plot(t_KF, -2*sig_KF(3,:), 'r--')
    xlabel("Time [sec]"); ylabel("Y [km]")
ax(4) = subplot(2,2,4);
    hold on; grid on;
    title("Ydot state error")
    plot(t_KF, stateError(4,:))
    plot(t_KF, 2*sig_KF(4,:), 'r--');
    plot(t_KF, -2*sig_KF(4,:), 'r--')
    xlabel("Time [sec]"); ylabel("Ydot [km/s]")
linkaxes(ax, 'x');

%% Measurement Perturbations
kTime = 1:length(filterData.dy_KF);

points = [55, 102, 149, 196, 244, 292, 340, 388, 436, 485, 533];

ax = [];

figure
sgtitle("Measurement Perturbations - Last Monte Carlo run")
ax(1) = subplot(3,1,1);
hold on; grid on;
title("\delta\rho vs. time")
for k = 1:length(filterData.dy_KF)
    if ~isempty(filterData.dy_KF{k})
        if k == 76
            a = plot(kTime(k), filterData.dy_KF{k}(1), 'b.');
            b = plot(kTime(k), filterData.yMeas_KF{k}(1), 'r.');
            c = plot(kTime(k), filterData.yNom_KF{k}(1), 'm.');
        else
            plot(kTime(k), filterData.dy_KF{k}(1), 'b.');
            plot(kTime(k), filterData.yMeas_KF{k}(1), 'r.')
            plot(kTime(k), filterData.yNom_KF{k}(1), 'm.')
        end
    end
end
plot(kTime(points), zeros(1,length(points)), 'k.', 'MarkerSize', 15);
legend([a,b,c], ["dy", "yMeas", "yNom"], 'location', 'bestoutside')

ax(2) = subplot(3,1,2);
hold on; grid on;
title("\delta\rho_{Dot} vs. time")
for k = 1:length(filterData.dy_KF)
    if ~isempty(filterData.dy_KF{k})
        if k == 76
            a = plot(kTime(k), filterData.dy_KF{k}(2), 'b.');
            b = plot(kTime(k), filterData.yMeas_KF{k}(2), 'r.');
            c = plot(kTime(k), filterData.yNom_KF{k}(2), 'm.');
        else
            plot(kTime(k), filterData.dy_KF{k}(2), 'b.');
            plot(kTime(k), filterData.yMeas_KF{k}(2), 'r.')
            plot(kTime(k), filterData.yNom_KF{k}(2), 'm.')
        end
    end
end
plot(kTime(points), zeros(1,length(points)), 'k.', 'MarkerSize', 15);
legend([a,b,c], ["dy", "yMeas", "yNom"], 'location', 'bestoutside')

ax(3) = subplot(3,1,3);
hold on; grid on;
title("\delta\phi vs. time")
for k = 1:length(filterData.dy_KF)
    if ~isempty(filterData.dy_KF{k})
        if k == 76
            a = plot(kTime(k), filterData.dy_KF{k}(3), 'b.');
            b = plot(kTime(k), filterData.yMeas_KF{k}(3), 'r.');
            c = plot(kTime(k), filterData.yNom_KF{k}(3), 'm.');
        else
            plot(kTime(k), filterData.dy_KF{k}(3), 'b.');
            plot(kTime(k), filterData.yMeas_KF{k}(3), 'r.')
            plot(kTime(k), filterData.yNom_KF{k}(3), 'm.')
        end
    end
end
plot(kTime(points), zeros(1,length(points)), 'k.', 'MarkerSize', 15);
legend([a,b,c], ["dy", "yMeas", "yNom"], 'location', 'bestoutside')

linkaxes(ax, 'x')


%% Implement on given data

filterData = LKF(ydata, filterParams, filterInit);
t_KF = filterData.t_KF;
x_KF = filterData.x_KF;
sig_KF = filterData.sig_KF;

figure
sgtitle("LKF Recovered States from Given Data")
ax(1) = subplot(2,2,1);
    hold on; grid on;
    title("X state")
    plot(t_KF, x_KF(1,:))
    % plot(t_KF, nomSolution(1,:), 'r--')
    xlabel("Time [sec]"); ylabel("X [km]")
ax(2) = subplot(2,2,2);
    hold on; grid on;
    title("Xdot state")
    plot(t_KF, x_KF(2,:))
    % plot(t_KF, nomSolution(2,:), 'r--')
    xlabel("Time [sec]"); ylabel("Xdot [km/s]")
ax(3) = subplot(2,2,3);
    hold on; grid on;
    title("Y state")
    plot(t_KF, x_KF(3,:))
    % plot(t_KF, nomSolution(3,:), 'r--')
    xlabel("Time [sec]"); ylabel("Y [km]")
ax(4) = subplot(2,2,4);
    hold on; grid on;
    title("Ydot state")
    plot(t_KF, x_KF(4,:))
    % plot(t_KF, nomSolution(4,:), 'r--')
    xlabel("Time [sec]"); ylabel("Ydot [km/s]")
linkaxes(ax, 'x');

figure
sgtitle("LKF Recovered 2\sigma bounds from Given Data")
ax(1) = subplot(2,2,1);
    hold on; grid on;
    title("X state bounds")
    plot(t_KF, 2*sig_KF(1,:), 'r--')
    plot(t_KF, -2*sig_KF(1,:), 'r--')
    xlabel("Time [sec]"); ylabel("X 2\sigma [km]")
ax(2) = subplot(2,2,2);
    hold on; grid on;
    title("Xdot state bounds")
    plot(t_KF, 2*sig_KF(2,:), 'r--')
    plot(t_KF, -2*sig_KF(2,:), 'r--')
    xlabel("Time [sec]"); ylabel("Xdot 2\sigma [km/s]")
ax(3) = subplot(2,2,3);
    hold on; grid on;
    title("Y state bounds")
    plot(t_KF, 2*sig_KF(3,:), 'r--')
    plot(t_KF, -2*sig_KF(3,:), 'r--')
    xlabel("Time [sec]"); ylabel("Y 2\sigma[km]")
ax(4) = subplot(2,2,4);
    hold on; grid on;
    title("Ydot state bounds")
    plot(t_KF, 2*sig_KF(4,:), 'r--')
    plot(t_KF, -2*sig_KF(4,:), 'r--')
    xlabel("Time [sec]"); ylabel("Ydot 2\sigma [km/s]")
linkaxes(ax, 'x');



