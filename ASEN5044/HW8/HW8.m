%% ASEN 5044 HW 8 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Problem 1
    % Helper functions
A_func = @(Omega) [
                        0 1     0 0;
                        0 0     0 -Omega;
                        0 0     0 1;
                        0 Omega 0 0
                  ];

F_func = @(Omega, dt)   [
                            1   sin(Omega*dt)/Omega     0   -(1-cos(Omega*dt))/Omega;
                            0   cos(Omega*dt)           0   -sin(Omega*dt);
                            0   (1-cos(Omega*dt))/Omega 1   sin(Omega*dt)/Omega;
                            0   sin(Omega*dt)           0   cos(Omega*dt)
                        ];

    % Setup
q_w = 10; % (m/s)^2
dt = 0.5; % sec

W = q_w*[
          2 0.05; 
          0.05 0.5
        ];

    % Aircraft A
Gamma_a = [
            0 0;
            1 0;
            0 0;
            0 1
          ];
Omega_a = 0.045; % rad/s
A_a = A_func(Omega_a);
F_a = F_func(Omega_a, dt)

Z_a = dt*[
            -A_a                Gamma_a*W*Gamma_a';
            zeros(size(A_a))    A_a'
         ];
expMat_a = expm(Z_a);

Q_a = F_a*expMat_a(1:4, 5:8)

    % Aircraft B
Gamma_b = [
            0 0;
            1 0;
            0 0;
            0 1
          ];
Omega_b = -0.045; % rad/s
A_b = A_func(Omega_b);
F_b = F_func(Omega_b, dt)

Z_b = dt*[
            -A_b                Gamma_b*W*Gamma_b';
            zeros(size(A_b))    A_b'
         ];
expMat_b = expm(Z_b);

Q_b = F_b*expMat_b(1:4, 5:8)

%% Problem 2 setup
rng(100);

%% Problem 2a
load("hw8problemdata.mat")

    % Setup
H = [
        1 0 0 0;
        0 0 1 0
    ];

R_a = [
        20   0.05;
        0.05 20
      ];

    % Make noisy samples
mv = zeros(2,1); Sv = chol(R_a, 'lower'); q = mvnrnd(zeros(1,2), eye(2), size(xasingle_truth, 2))';
v = mv + Sv*q; % Create noise samples for k = 1:100/dt (100 seconds)

    % Recreate measurements
y = [];
for k = 2:size(xasingle_truth,2)
    y = [y, H*xasingle_truth(:,k) + v(:,k)];
end

% Plot first 20 seconds of measurements
tStop = 20; idx = 1:tStop/dt; time = idx*dt;
titleText = sprintf("y_A(k) vs. time, first %.0f seconds", tStop);
figure;
sgtitle(titleText)
subplot(2,1,1)
    hold on; grid on;
    plot(time, y(1,idx));
    xlabel("Time [sec]"); ylabel("\xi_A [m]")
subplot(2,1,2)
    hold on; grid on;
    plot(time, y(2,idx));
    xlabel("Time [sec]"); ylabel("\eta_A [m]")

%% Problem 2b
    % Setup
mu_a0 = [0; 85*cos(pi/4); 0; -85*sin(pi/4)];
P_a0 = 900*diag([10 2 10 2]);
F = F_a;
Q = Q_a;

    % Kalman Filter
t_KF = [];
x_KF = [];
P_KF = zeros(size(F,1), size(F,2), size(y,2));
innovation_KF = [];
Sk_KF = zeros(size(H,1), size(H,1), size(y,2));
sig_KF = [];
x_kPlus = mu_a0;
P_kPlus = P_a0;
for k = 1:size(y,2)
    % Matrix assignment
    H_kp1 = H;
    R_kp1 = R_a;

    % Time Update / Prediction Step
    x_kp1Minus = F*x_kPlus;
    P_kp1Minus = F*P_kPlus*F' + Q;
    K_kp1 = P_kp1Minus*H_kp1'*(H_kp1*P_kp1Minus*H_kp1' + R_kp1)^-1;

    % Measurement Update / Correction Step
    x_kp1Plus = x_kp1Minus + K_kp1*(y(:,k) - H_kp1*x_kp1Minus); % y starts at y_1, not y_0
    P_kp1Plus = (eye(size(F))-K_kp1*H_kp1)*P_kp1Minus;

    % Save outputs
    t_KF = [t_KF; k*dt];
    x_KF = [x_KF, x_kp1Plus];
    P_KF(:,:,k) = P_kp1Plus;
    innovation_KF = [innovation_KF, y(:,k) - H_kp1*x_kp1Minus];
    Sk_KF(:,:,k) = H_kp1*P_kp1Minus*H' + R_kp1;
    sig_KF = [sig_KF, [sqrt(P_KF(1,1,k)); sqrt(P_KF(2,2,k)); sqrt(P_KF(3,3,k)); sqrt(P_KF(4,4,k))]];

    % Update for next run
    x_kPlus = x_kp1Plus;
    P_kPlus = P_kp1Plus;

end

    % Plot state errors
figure; ax = [];
sgtitle("Kalman Filter State Errors vs. Time")
ax(1) = subplot(4,1,1);
    hold on; grid on;
    plot(t_KF, xasingle_truth(1,2:end) - x_KF(1,:), 'b-')
    plot(t_KF, 2*sig_KF(1,:), 'r--')
    plot(t_KF, -2*sig_KF(1,:), 'r--')
    xlabel("Time [sec]"); ylabel("\xi_A error [m]")
ax(2) = subplot(4,1,2);
    hold on; grid on;
    plot(t_KF, xasingle_truth(2,2:end) - x_KF(2,:), 'b-')
    plot(t_KF, 2*sig_KF(2,:), 'r--')
    plot(t_KF, -2*sig_KF(2,:), 'r--')
    xlabel("Time [sec]"); ylabel("\textbf{$\dot{\xi}_A$ error [m/s]}", 'Interpreter', 'latex')
ax(3) = subplot(4,1,3);
    hold on; grid on;
    plot(t_KF, xasingle_truth(3,2:end) - x_KF(3,:), 'b-')
    plot(t_KF, 2*sig_KF(3,:), 'r--')
    plot(t_KF, -2*sig_KF(3,:), 'r--')
    xlabel("Time [sec]"); ylabel("\eta_A error [m]")
ax(4) = subplot(4,1,4);
    hold on; grid on;
    plot(t_KF, xasingle_truth(4,2:end) - x_KF(4,:), 'b-')
    plot(t_KF, 2*sig_KF(4,:), 'r--')
    plot(t_KF, -2*sig_KF(4,:), 'r--')
    xlabel("Time [sec]"); ylabel("\textbf{$\dot{\eta}_A$ error [m/s]}", 'Interpreter', 'latex')
linkaxes(ax, 'x')

%     % Plot state estimates
% figure
% sgtitle("Kalman Filter State Estimates vs. Time")
% ax(1) = subplot(4,1,1);
%     hold on; grid on;
%     plot(t_KF, x_KF(1,:), 'b-')
%     plot(t_KF, xasingle_truth(1,2:end), 'b--')
%     plot(t_KF, x_KF(1,:) + 2*sig_KF(1,:), 'r--')
%     plot(t_KF, x_KF(1,:) - 2*sig_KF(1,:), 'r--')
%     xlabel("Time [sec]"); ylabel("\xi_A [m]")
% ax(2) = subplot(4,1,2);
%     hold on; grid on;
%     plot(t_KF, x_KF(2,:), 'b-')
%     plot(t_KF, xasingle_truth(2,2:end), 'b--')
%     plot(t_KF, x_KF(2,:) + 2*sig_KF(2,:), 'r--')
%     plot(t_KF, x_KF(2,:) - 2*sig_KF(2,:), 'r--')
%     xlabel("Time [sec]"); ylabel("\xiDot_A [m]")
% ax(3) = subplot(4,1,3);
%     hold on; grid on;
%     plot(t_KF, x_KF(3,:), 'b-')
%     plot(t_KF, xasingle_truth(3,2:end), 'b--')
%     plot(t_KF, x_KF(3,:) + 2*sig_KF(3,:), 'r--')
%     plot(t_KF, x_KF(3,:) - 2*sig_KF(3,:), 'r--')
%     xlabel("Time [sec]"); ylabel("\eta_A [m]")
% ax(4) = subplot(4,1,4);
%     hold on; grid on;
%     plot(t_KF, x_KF(4,:), 'b-')
%     plot(t_KF, xasingle_truth(4,2:end), 'b--')
%     plot(t_KF, x_KF(4,:) + 2*sig_KF(4,:), 'r--')
%     plot(t_KF, x_KF(4,:) - 2*sig_KF(4,:), 'r--')
%     xlabel("Time [sec]"); ylabel("\etaDot_A [m]")
% linkaxes(ax, 'x')

%% Problem 2 NEES/NIS tests
alpha = 0.05;
N = 1+0*size(x_KF, 2);
n = size(F,1);
p = size(H,1);
r1_NEES = chi2inv(alpha/2, N*n)./N;
r2_NEES = chi2inv(1-(alpha/2), N*n)./N;
r1_NIS = chi2inv(alpha/2, N*p)./N;
r2_NIS = chi2inv(1-(alpha/2), N*p)./N;

epsilon_x = [];
epsilon_y = [];
for k = 1:size(x_KF, 2)
    stateError = xasingle_truth(:,k+1) - x_KF(:,k);
    epsilon_x = [epsilon_x; stateError'*(P_KF(:,:,k)^-1)*stateError];
    epsilon_y = [epsilon_y; innovation_KF(:,k)'*(Sk_KF(:,:,k)^-1)*innovation_KF(:,k)];
end

figure
title("NEES Statistic (N = 1, \alpha = 0.05)")
hold on; grid on;
plot(epsilon_x, 'ro')
yline(r1_NEES, 'r--')
yline(r2_NEES, 'r--')
xlabel("k"); ylabel("$\bar{\epsilon}_{x,k}$", 'Interpreter', 'latex')
legend("NEES at time k", "r_1 bound", "r_2 bound")

figure
title("NIS Statistic (N = 1, \alpha = 0.05)")
hold on; grid on;
plot(epsilon_y, 'ro')
yline(r1_NIS, 'r--')
yline(r2_NIS, 'r--')
xlabel("k"); ylabel("$\bar{\epsilon}_{y,k}$", 'Interpreter', 'latex')

%% Test function
filter2b = KalmanFilter(y, F_a, Q_a, H, R_a, mu_a0, P_a0, dt, zeros(2,1));

%% Problem 3a
    % Setup
R_d = [
            10   0.15;
            0.15 10
      ];

    % Initial conditions
mu_a0 = [0; 85*cos(pi/4); 0; -85*sin(pi/4)];
P_a0 = 900*diag([10, 2, 10, 2]);

mu_b0 = [3200; 85*cos(pi/4); 3200; -85*sin(pi/4)];
P_b0 = 900*diag([11, 4, 11, 4]);

    % New matrices
F_s = blkdiag(F_a, F_b);

A_s = blkdiag(A_a, A_b);
Gamma_s = blkdiag(Gamma_a, Gamma_b);
W_s = blkdiag(W, W);
Z_s = dt*[
            -A_s                Gamma_s*W_s*Gamma_s';
            zeros(size(A_s))    A_s'
         ];
matExp_s = expm(Z_s);
Q_s = F_s*matExp_s(1:8, 9:16);

H_s = [
            H   zeros(size(H));
            H   -H
      ];

R_s = blkdiag(R_a, R_d);

x_s0 = [mu_a0; mu_b0];

P_s0 = blkdiag(P_a0, P_b0);

    % Simulate noisy measurements
mv_a = zeros(2,1); Sv_a = chol(R_a, 'lower'); 
q_a = mvnrnd(zeros(1,2), eye(2), size(xadouble_truth, 2))';
v_a = mv_a + Sv_a*q_a; % Create noise samples for k = 1:100/dt (100 seconds)

mv_d = zeros(2,1); Sv_d = chol(R_d, 'lower');
q_d = mvnrnd(zeros(1,2), eye(2), size(xadouble_truth, 2))';
v_d = mv_d + Sv_d*q_d; % Create noise samples for k = 1:100/dt (100 seconds)

v_s = [v_a; v_d];

y_s = [];
for k = 2:size(xadouble_truth, 2)
    y_s = [y_s, H_s*[xadouble_truth(:,k); xbdouble_truth(:,k)] + v_s(:,k)];
end

    % Run Kalman Filter
filter3a = KalmanFilter(y_s, F_s, Q_s, H_s, R_s, x_s0, P_s0, dt, zeros(2,1));

    % Extract results
t_KF = filter3a.t_KF;
x_KF = filter3a.x_KF;
sig_KF = filter3a.sig_KF;

    % Compute state errors
stateError = [xadouble_truth(:,2:end); xbdouble_truth(:,2:end)] - x_KF;

    % Plot position errors with 2 sigma bounds
figure; ax = [];
t = tiledlayout(2,1);
title(t, "Position Errors for Aircraft A - Ground Station and Transponder")
ax(1) = nexttile;
    hold on; grid on;
    title("\xi_A error vs. time")
    error = plot(t_KF, stateError(1,:), 'b-');
    sigma = plot(t_KF, 2*sig_KF(1,:), 'r--');
    plot(t_KF, -2*sig_KF(1,:), 'r--')
    xlabel("Time [sec]"); ylabel("\xi_A error [m]")
ax(2) = nexttile;
    hold on; grid on;
    title("\eta_A vs. time")
    plot(t_KF, stateError(3,:), 'b-')
    plot(t_KF, 2*sig_KF(3,:), 'r--')
    plot(t_KF, -2*sig_KF(3,:), 'r--')
    xlabel("Time [sec]"); ylabel("\eta_A error [m]")
linkaxes(ax, 'x')
legend([error, sigma], ["State error", "2\sigma bounds"], 'location', 'bestoutside')

figure; ax = [];
t = tiledlayout(2,1);
title(t, "Position Errors for Aircraft B - Ground Station and Transponder")
ax(1) = nexttile;
    hold on; grid on;
    title("\xi_B error vs. time")
    error = plot(t_KF, stateError(5,:), 'b-');
    sigma = plot(t_KF, 2*sig_KF(5,:), 'r--');
    plot(t_KF, -2*sig_KF(5,:), 'r--')
    xlabel("Time [sec]"); ylabel("\xi_B error [m]")
ax(2) = nexttile;
    hold on; grid on;
    title("\eta_B vs. time")
    plot(t_KF, stateError(7,:), 'b-')
    plot(t_KF, 2*sig_KF(7,:), 'r--')
    plot(t_KF, -2*sig_KF(7,:), 'r--')
    xlabel("Time [sec]"); ylabel("\eta_B error [m]")
linkaxes(ax, 'x')
legend([error, sigma], ["State error", "2\sigma bounds"], 'location', 'bestoutside')

%     % Plot position estimates with 2 sigma bounds
% ax = [];
% figure
% sgtitle("Position Estimates for Aircraft A")
% ax(1) = subplot(2,1,1);
%     hold on; grid on;
%     title("\xi_A vs. time")
%     state = plot(t_KF, x_KF(1,:), 'b-');
%     data = plot(t_KF, xadouble_truth(1,2:end), 'b--');
%     sigma = plot(t_KF, x_KF(1,:) + 2*sig_KF(1,:), 'r--');
%     plot(t_KF, x_KF(1,:) - 2*sig_KF(1,:), 'r--')
%     xlabel("Time [sec]"); ylabel("\xi_A [m]")
% ax(2) = subplot(2,1,2);
%     hold on; grid on;
%     title("\eta_A vs. time")
%     plot(t_KF, x_KF(3,:), 'b-')
%     plot(t_KF, xadouble_truth(3,2:end), 'b--');
%     plot(t_KF, x_KF(3,:) + 2*sig_KF(3,:), 'r--')
%     plot(t_KF, x_KF(3,:) - 2*sig_KF(3,:), 'r--')
%     xlabel("Time [sec]"); ylabel("\eta_A [m]")
% linkaxes(ax, 'x')
% legend([state, data, sigma], ["Estimated state", "Measured data", "2\sigma bounds"], 'location', 'best')
% 
% ax = [];
% figure
% sgtitle("Position Estimates for Aircraft B")
% ax(1) = subplot(2,1,1);
%     hold on; grid on;
%     title("\xi_B vs. time")
%     state = plot(t_KF, x_KF(5,:), 'b-');
%     data = plot(t_KF, xbdouble_truth(1,2:end), 'b--');
%     sigma = plot(t_KF, x_KF(5,:) + 2*sig_KF(5,:), 'r--');
%     plot(t_KF, x_KF(5,:) - 2*sig_KF(5,:), 'r--')
%     xlabel("Time [sec]"); ylabel("\xi_B [m]")
% ax(2) = subplot(2,1,2);
%     hold on; grid on;
%     title("\eta_B vs. time")
%     plot(t_KF, x_KF(7,:), 'b-')
%     plot(t_KF, xbdouble_truth(3,2:end), 'b--');
%     plot(t_KF, x_KF(7,:) + 2*sig_KF(7,:), 'r--')
%     plot(t_KF, x_KF(7,:) - 2*sig_KF(7,:), 'r--')
%     xlabel("Time [sec]"); ylabel("\eta_A [m]")
% linkaxes(ax, 'x')
% legend([state, data, sigma], ["Estimated state", "Measured data", "2\sigma bounds"], 'location', 'best')

%% Problem 3b
    % New matrices
H_s = [H, -H];
R_s = R_d;

    % New measurements
v_s = v_d; % Same rng, no need to regenerate

y_s = [];
for k = 2:size(xadouble_truth, 2)
    y_s = [y_s, H_s*[xadouble_truth(:,k); xbdouble_truth(:,k)] + v_s(:,k)];
end

    % Kalman Filter
filter3b = KalmanFilter(y_s, F_s, Q_s, H_s, R_s, x_s0, P_s0, dt, zeros(2,1));

    % Extract results
t_KF = filter3b.t_KF;
x_KF = filter3b.x_KF;
sig_KF = filter3b.sig_KF;

    % Compute state errors
stateError = [xadouble_truth(:,2:end); xbdouble_truth(:,2:end)] - x_KF;

    % Plot position errors with 2 sigma bounds
figure; ax = [];
t = tiledlayout(2,1);
title(t, "Position Errors for Aircraft A - Transponder Only")
ax(1) = nexttile;
    hold on; grid on;
    title("\xi_A error vs. time")
    error = plot(t_KF, stateError(1,:), 'b-');
    sigma = plot(t_KF, 2*sig_KF(1,:), 'r--');
    plot(t_KF, -2*sig_KF(1,:), 'r--')
    xlabel("Time [sec]"); ylabel("\xi_A error [m]")
ax(2) = nexttile;
    hold on; grid on;
    title("\eta_A vs. time")
    plot(t_KF, stateError(3,:), 'b-')
    plot(t_KF, 2*sig_KF(3,:), 'r--')
    plot(t_KF, -2*sig_KF(3,:), 'r--')
    xlabel("Time [sec]"); ylabel("\eta_A error [m]")
linkaxes(ax, 'x')
legend([error, sigma], ["State error", "2\sigma bounds"], 'location', 'bestoutside')

figure; ax = [];
t = tiledlayout(2,1);
title(t,"Position Errors for Aircraft B - Transponder Only")
ax(1) = nexttile;
    hold on; grid on;
    title("\xi_B error vs. time")
    error = plot(t_KF, stateError(5,:), 'b-');
    sigma = plot(t_KF, 2*sig_KF(5,:), 'r--');
    plot(t_KF, -2*sig_KF(5,:), 'r--')
    xlabel("Time [sec]"); ylabel("\xi_B error [m]")
ax(2) = nexttile;
    hold on; grid on;
    title("\eta_B vs. time")
    plot(t_KF, stateError(7,:), 'b-')
    plot(t_KF, 2*sig_KF(7,:), 'r--')
    plot(t_KF, -2*sig_KF(7,:), 'r--')
    xlabel("Time [sec]"); ylabel("\eta_B error [m]")
linkaxes(ax, 'x')
legend([error, sigma], ["State error", "2\sigma bounds"], 'location', 'bestoutside')

%% Problem 3c
    % Measurement steps shouldn't matter, are ignored in this problem
    % New matrices
H_s = [
            H   zeros(size(H));
            H   -H
      ];

R_s = blkdiag(R_a, R_d);

    % New measurements
v_s = [v_a; v_d];

y_s = [];
for k = 2:size(xadouble_truth, 2)
    y_s = [y_s, H_s*[xadouble_truth(:,k); xbdouble_truth(:,k)] + v_s(:,k)];
end

    % Set debug - want only prediction in this problem
noMeas = true; noPred = false; debug = [noMeas; noPred];

    % Kalman Filter
filter3c = KalmanFilter(y_s, F_s, Q_s, H_s, R_s, x_s0, P_s0, dt, debug);

    % Extract results
t_KF = filter3c.t_KF;
x_KF = filter3c.x_KF;
sig_KF = filter3c.sig_KF;

    % Compute state errors
stateError = [xadouble_truth(:,2:end); xbdouble_truth(:,2:end)] - x_KF;

    % Plot position errors with 2 sigma bounds
figure; ax = [];
t = tiledlayout(2,1);
title(t, "Position Errors for Aircraft A - Prediction only")
ax(1) = nexttile;
    hold on; grid on;
    title("\xi_A error vs. time")
    error = plot(t_KF, stateError(1,:), 'b-');
    sigma = plot(t_KF, 2*sig_KF(1,:), 'r--');
    plot(t_KF, -2*sig_KF(1,:), 'r--')
    xlabel("Time [sec]"); ylabel("\xi_A error [m]")
ax(2) = nexttile;
    hold on; grid on;
    title("\eta_A vs. time")
    plot(t_KF, stateError(3,:), 'b-')
    plot(t_KF, 2*sig_KF(3,:), 'r--')
    plot(t_KF, -2*sig_KF(3,:), 'r--')
    xlabel("Time [sec]"); ylabel("\eta_A error [m]")
linkaxes(ax, 'x')
legend([error, sigma], ["State error", "2\sigma bounds"], 'location', 'bestoutside')

figure; ax = [];
t = tiledlayout(2,1);
title(t,"Position Errors for Aircraft B - Prediction Only")
ax(1) = nexttile;
    hold on; grid on;
    title("\xi_B error vs. time")
    error = plot(t_KF, stateError(5,:), 'b-');
    sigma = plot(t_KF, 2*sig_KF(5,:), 'r--');
    plot(t_KF, -2*sig_KF(5,:), 'r--')
    xlabel("Time [sec]"); ylabel("\xi_B error [m]")
ax(2) = nexttile;
    hold on; grid on;
    title("\eta_B vs. time")
    plot(t_KF, stateError(7,:), 'b-')
    plot(t_KF, 2*sig_KF(7,:), 'r--')
    plot(t_KF, -2*sig_KF(7,:), 'r--')
    xlabel("Time [sec]"); ylabel("\eta_B error [m]")
linkaxes(ax, 'x')
legend([error, sigma], ["State error", "2\sigma bounds"], 'location', 'bestoutside')
