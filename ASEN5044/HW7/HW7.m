%% ASEN 5044 HW 7 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Problem 2
% Constants
dt = 0.5; % sec
Omega = 0.045; % rad/s
mu0 = [0; 85*cos(pi/4); 0; -85*sin(pi/4)];
Pa0 = diag([10 2 10 2]);

% Define DT system
A = [0 1 0 0; 0 0 0 -Omega; 0 0 0 1; 0 Omega 0 0];
matExp = expm(A*dt); % Check F

F_func = @(Omega, dt)   [
                            1   sin(Omega*dt)/Omega     0   -(1-cos(Omega*dt))/Omega;
                            0   cos(Omega*dt)           0   -sin(Omega*dt);
                            0   (1-cos(Omega*dt))/Omega 1   sin(Omega*dt)/Omega;
                            0   sin(Omega*dt)           0   cos(Omega*dt)
                        ];

% Propagate mean and covariance matrices
t = 1:300;
F = F_func(Omega, dt);
Mu = zeros(size(F,1),1,length(t));
P = zeros(size(F,1),size(F,2),length(t));
twoSigma = zeros(size(F,1),1,length(t));
for k = t
    kMath = k - 1;
    Mu(:,1,k) = (F^kMath)*mu0;
    P(:,:,k) = (F^kMath)*Pa0*(F^kMath)';
    twoSigma(:,:,k) = 2*[sqrt(P(1,1,k)); round(sqrt(P(2,2,k)),6); sqrt(P(3,3,k)); round(sqrt(P(4,4,k)),6)];
end

figure;
sgtitle("\mu(k) vs. time")
ax(1) = subplot(4,1,1);
    hold on; grid on;
    title("\xi vs. time")
    muXi = plot(t, reshape(Mu(1,1,:), 1, []));
    sigmaXi = plot(t, reshape(Mu(1,1,:) + twoSigma(1,1,:), 1 ,[]), 'r--');
    plot(t, reshape(Mu(1,1,:) - twoSigma(1,1,:), 1, []), 'r--')
    xlabel("Time [k]"); ylabel("\xi [m]")
    legend([muXi, sigmaXi], ["\mu(k)", "+/- 2 sigma bound"], 'location', 'best')
ax(2) = subplot(4,1,2);
    hold on; grid on;
    title("\xiDot vs. time")
    plot(t, reshape(Mu(2,1,:), 1, []));
    plot(t, reshape(Mu(2,1,:) + twoSigma(2,1,:), 1 ,[]), 'r--');
    plot(t, reshape(Mu(2,1,:) - twoSigma(2,1,:), 1, []), 'r--')
    xlabel("Time [k]"); ylabel("\xiDot [m/s]")
ax(3) = subplot(4,1,3);
    hold on; grid on;
    title("\eta vs. time")
    plot(t, reshape(Mu(3,1,:), 1, []));
    plot(t, reshape(Mu(3,1,:) + twoSigma(3,1,:), 1 ,[]), 'r--');
    plot(t, reshape(Mu(3,1,:) - twoSigma(3,1,:), 1, []), 'r--')
    xlabel("Time [k]"); ylabel("\eta [m]")
ax(4) = subplot(4,1,4);
    hold on; grid on;
    title("\etaDot vs. time")
    plot(t, reshape(Mu(4,1,:), 1, []));
    plot(t, reshape(Mu(4,1,:) + twoSigma(4,1,:), 1 ,[]), 'r--');
    plot(t, reshape(Mu(4,1,:) - twoSigma(4,1,:), 1, []), 'r--')
    xlabel("Time [k]"); ylabel("\etaDot [m/s]")

figure;
sgtitle("2\sigma(k) vs. time")
ax(5) = subplot(4,1,1);
    hold on; grid on;
    title("2\sigma_\xi vs. time")
    plot(t, reshape(twoSigma(1,1,:), 1, []), 'r--')
    xlabel("Time [k]"); ylabel("2\sigma_\xi [m]")
ax(6) = subplot(4,1,2);
    hold on; grid on;
    title("2\sigma_{\xiDot} vs. time")
    plot(t, reshape(twoSigma(2,1,:), 1, []), 'r--')
    xlabel("Time [k]"); ylabel("2\sigma_{\xiDot} [m/s]")
ax(7) = subplot(4,1,3);
    hold on; grid on;
    title("2\sigma_\eta vs. time")
    plot(t, reshape(twoSigma(3,1,:), 1, []), 'r--')
    xlabel("Time [k]"); ylabel("2\sigma_\eta [m]")
ax(8) = subplot(4,1,4);
    hold on; grid on;
    title("2\sigma_{\etaDot} vs. time")
    plot(t, reshape(twoSigma(4,1,:), 1, []), 'r--')
    xlabel("Time [k]"); ylabel("2\sigma_{\etaDot} [m/s]")
    
linkaxes(ax, 'x')

    % Sanity check
if false
    figure
    hold on; grid on;
    title("Trajectory of \mu(k)")
    plot(reshape(Mu(1,1,:), 1, []), reshape(Mu(3,1,:), 1, []));
    startTraj = plot(Mu(1,1,1), Mu(3,1,1), 'g.', 'MarkerSize', 15);
    endTraj = plot(Mu(1,1,end), Mu(3,1,end), 'r.', 'MarkerSize', 15);
    xlabel("\xi [m]"); ylabel("\eta [m]");
    legend([startTraj, endTraj], ["Start of trajectory", "End of trajectory"], 'Location', 'best')
end

%% Problem 3
    % Problem params
dt = 0.5; % sec
t = 0:dt:150;
xiR = 100; % m
etaR = 100; % m

    % Aircraft a parameters
Omega_a = 0.045; % rad/s
mu_a0 = [0; 85*cos(pi/4); 0; -85*sin(pi/4)];
P_a0 = diag([10, 4, 10, 4]);
F_a = F_func(Omega_a, dt);

    % Aircraft b parameters
Omega_b = -0.045; % rad/s
mu_b0 = [3200; 85*cos(pi/4); 3200; -85*sin(pi/4)];
P_b0 = diag([11, 3.5, 11, 3.5]);
F_b = F_func(Omega_b, dt);

    % Collision probability parameters
M = [
        1 0 0 0;
        0 0 1 0
    ];
mu_rc = @(k) M*((F_a^k)*mu_a0 - (F_b^k)*mu_b0);
P_rc = @(k) M*((F_a^k)*P_a0*(F_a^k)' + (F_b^k)*P_b0*(F_b^k)')*M';

    % Simulate probability
p = zeros(length(t), 1);
for k = 1:length(t)
    kMath = k - 1;
    
    mu = mu_rc(kMath);
    P = P_rc(kMath);

    p(k) = mvncdf([-xiR; -etaR], [xiR; etaR], mu, P);
end

    % Simulate expected trajectories
mu_a = zeros(size(F_a,1), 1, length(t));
mu_b = zeros(size(F_b,1), 1, length(t));
for k = 1:length(t)
    kMath = k - 1;
    mu_a(:,1,k) = (F_a^k)*mu_a0;
    mu_b(:,1,k) = (F_b^k)*mu_b0;
end

figure
hold on; grid on;
title("Probability of collision vs. time")
plot(t, p)
xlabel("Time [sec]"); ylabel("Probability of collision")

[colP, colIdx] = maxk(p,8);

for k = 1:length(colIdx)
    fprintf("%.3f%% chance of collision at t = %.1f sec (k = %.0f)!!\n", 100*colP(k), dt*colIdx(k), colIdx(k))
end

%% Problem 3 Animation
figure
for k = 1:length(mu_a)
    clf; hold on; grid on;

    titleText = sprintf("Expected trajectories of aircraft a and b at t = %.1f sec", (k-1)*dt);
    title(titleText)

    plot(reshape(mu_a(1,1,:), 1, []), reshape(mu_a(3,1,:), 1, []), 'b-');
    plot(reshape(mu_b(1,1,:), 1, []), reshape(mu_b(3,1,:), 1, []), 'r-');
    startTraj_a = plot(mu_a(1,1,1), mu_a(3,1,1), 'go', 'MarkerSize', 15);
    endTraj_a = plot(mu_a(1,1,end), mu_a(3,1,end), 'ro', 'MarkerSize', 15);
    startTraj_b = plot(mu_b(1,1,1), mu_b(3,1,1), 'g*', 'MarkerSize', 15);
    endTraj_b = plot(mu_b(1,1,end), mu_b(3,1,end), 'r*', 'MarkerSize', 15);
    col_1 = plot(mu_a(1,1,colIdx(1)), mu_a(3,1,colIdx(1)), 'mx', 'Markersize', 15);
    col_2 = plot(mu_a(1,1,colIdx(2)), mu_a(3,1,colIdx(2)), 'cx', 'Markersize', 15);

    plot(reshape(mu_a(1,1,k), 1, []), reshape(mu_a(3,1,k), 1, []), 'k.', 'MarkerSize', 20);
    plot(reshape(mu_b(1,1,k), 1, []), reshape(mu_b(3,1,k), 1, []), 'k.', 'MarkerSize', 20);

    col_1Text = sprintf("%.3f%% chance of collision", 100*colP(1));
    col_2Text = sprintf("%.3f%% chance of collision", 100*colP(2));
    xlabel("\xi [m]"); ylabel("\eta [m]");
    legend([startTraj_a, endTraj_a, startTraj_b, endTraj_b, col_1, col_2], ["Start of aircraft a trajectory", "End of aircraft a trajectory", ...
            "Start of aircraft b trajectory", "End of aircraft b trajectory", col_1Text, col_2Text], 'Location', 'bestoutside')

    drawnow
end

