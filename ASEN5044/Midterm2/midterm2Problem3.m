%% ASEN 5044 Midterm 2 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
dt = 0.5; % sec

F_func = @(Omega, dt)   [
                            1   sin(Omega*dt)/Omega     0   -(1-cos(Omega*dt))/Omega;
                            0   cos(Omega*dt)           0   -sin(Omega*dt);
                            0   (1-cos(Omega*dt))/Omega 1   sin(Omega*dt)/Omega;
                            0   sin(Omega*dt)           0   cos(Omega*dt)
                        ];

%% Part b
Omega_a = 0.045; % rad/s
F_a = F_func(Omega_a, dt);

R_func = @(k)   [
                    75 + 12.5*sin(k/10)     7.5 + 25.5*sin(k/10)
                    7.5 + 25.5*sin(k/10)    75 + 12.5*cos(k/10)
                ];

H = [
        1 0 0 0;
        0 0 1 0
    ];

load("midterm2_problem3b.mat")

Rb = [];
HF = [];
y = [];
for k = 1:size(yaHist,2)
    Rb = blkdiag(Rb, R_func(k));
    HF = [HF; H*F_a^k];
    y = [y; yaHist(:,k)];
end

xHat0_a = ((HF'*(Rb^-1)*HF)^-1)*HF'*(Rb^-1)*y

Pls0 = (HF'*(Rb^-1)*HF)^-1

%% Part c
Omega_b = -0.045; % rad/s
F_b = F_func(Omega_b, dt);

RD = [
        8000 500;
        500  8000
     ];

load("midterm2_problem3c.mat")

xEst = [];
Pest = zeros(8,8,size(yaugHist,2));
twoSig_a0 = [];
twoSig_b0 = [];
xLS = zeros(8,1);
Pls = 9999*eye(8);
for k = 1:size(yaugHist,2)
    % Calculate time varying quantities
    Rk = blkdiag(R_func(k), RD);
    Hk = [
            H*F_a^k zeros(2,4);
            H*F_a^k -H*F_b^k
         ];

    % Propagate estimator
    Kk = Pls*Hk'*(Rk + Hk*Pls*Hk')^-1;
    xLS = xLS + Kk*(yaugHist(:,k)-Hk*xLS);
    Pls = (eye(8) - Kk*Hk)*Pls*(eye(8) - Kk*Hk)' + Kk*Rk*Kk';

    % Save estimates
    xEst = [xEst, xLS];
    Pest(:,:,k) = Pls;
    twoSig_a0 = [twoSig_a0, 2*[sqrt(Pest(1,1,k)); sqrt(Pest(2,2,k)); sqrt(Pest(3,3,k)); sqrt(Pest(4,4,k))]];
    twoSig_b0 = [twoSig_b0, 2*[sqrt(Pest(5,5,k)); sqrt(Pest(6,6,k)); sqrt(Pest(7,7,k)); sqrt(Pest(8,8,k))]];

end

x_a0 = xEst(1:4, end)
x_b0 = xEst(5:8, end)

figure
sgtitle("RLLS $\hat{x}_A(0)$ vs. k", 'Interpreter', 'latex')
subplot(4,1,1)
    hold on; grid on;
    title("\xi(0) vs. k")
    plot(xEst(1,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\xi(0) [m]")
subplot(4,1,2)
    hold on; grid on;
    title("\xiDot(0) vs. k")
    plot(xEst(2,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\xiDot(0) [m/s]")
subplot(4,1,3)
    hold on; grid on;
    title("\eta(0) vs. k")
    plot(xEst(3,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\eta(0) [m]")
subplot(4,1,4)
    hold on; grid on;
    title("\etaDot(0) vs. k")
    plot(xEst(4,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\etaDot(0) [m/s]")

figure
sgtitle("RLLS $\hat{x}_A(0)$ upper $2\sigma$ bound vs. k", 'Interpreter', 'latex')
subplot(4,1,1)
    hold on; grid on;
    title("2\sigma_{\xi(0)} vs. k")
    plot(twoSig_a0(1,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_{\xi(0)} [m]")
subplot(4,1,2)
    hold on; grid on;
    title("2\sigma_{\xiDot(0)} vs. k")
    plot(twoSig_a0(2,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_{\xiDot(0)} [m/s]")
subplot(4,1,3)
    hold on; grid on;
    title("2\sigma_{\eta(0)} vs. k")
    plot(twoSig_a0(3,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_{\eta(0)} [m]")
subplot(4,1,4)
    hold on; grid on;
    title("2\sigma_{\etaDot(0)} vs. k")
    plot(twoSig_a0(4,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_{\etaDot(0)} [m/s]")

figure
sgtitle("RLLS $\hat{x}_B(0)$ vs. k", 'Interpreter', 'latex')
subplot(4,1,1)
    hold on; grid on;
    title("\xi(0) vs. k")
    plot(xEst(5,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\xi(0) [m]")
subplot(4,1,2)
    hold on; grid on;
    title("\xiDot(0) vs. k")
    plot(xEst(6,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\xiDot(0) [m/s]")
subplot(4,1,3)
    hold on; grid on;
    title("\eta(0) vs. k")
    plot(xEst(7,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\eta(0) [m]")
subplot(4,1,4)
    hold on; grid on;
    title("\etaDot(0) vs. k")
    plot(xEst(8,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\etaDot(0) [m/s]")

figure
sgtitle("RLLS $\hat{x}_B(0)$ upper $2\sigma$ bound vs. k", 'Interpreter', 'latex')
subplot(4,1,1)
    hold on; grid on;
    title("2\sigma_{\xi(0)} vs. k")
    plot(twoSig_b0(1,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_{\xi(0)} [m]")
subplot(4,1,2)
    hold on; grid on;
    title("2\sigma_{\xiDot(0)} vs. k")
    plot(twoSig_b0(2,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_{\xiDot(0)} [m/s]")
subplot(4,1,3)
    hold on; grid on;
    title("2\sigma_{\eta(0)} vs. k")
    plot(twoSig_b0(3,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_{\eta(0)} [m]")
subplot(4,1,4)
    hold on; grid on;
    title("2\sigma_{\etaDot(0)} vs. k")
    plot(twoSig_b0(4,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_{\etaDot(0)} [m/s]")

%% Part d

% Estimated aircraft states and uncertainties
x_a = [];
twoSig_a = [];
P_a0 = Pest(1:4, 1:4, end);
x_b = [];
twoSig_b = [];
P_b0 = Pest(5:8, 5:8, end);
x_bData = [];
for k = 1:1*size(yaugHist,2)
    x_a = [x_a, F_a^k*x_a0];
    P_a = F_a^k*P_a0*(F_a^k)';
    twoSig_a = [twoSig_a, 2*[sqrt(P_a(1,1)); sqrt(P_a(2,2)); sqrt(P_a(3,3)); sqrt(P_a(4,4))]];
    
    x_b = [x_b, F_b^k*x_b0];
    P_b = F_b^k*P_b0*(F_b^k)';
    twoSig_b = [twoSig_b, 2*[sqrt(P_b(1,1)); sqrt(P_b(2,2)); sqrt(P_b(3,3)); sqrt(P_b(4,4))]];

    if k <= size(yaugHist,2) % Allow for dynamic propagation
        x_bData = [x_bData, yaugHist(1:2,k) - yaugHist(3:4,k)]; % y_D = r_A - r_B + v_D -> r_B = r_A - y_D + v_D
    end
end

% Plot aircraft A's state
figure
sgtitle("$\hat{x}_A(k)$ vs. k", 'Interpreter', 'latex')
subplot(4,1,1)
    hold on; grid on;
    title("\xi vs. k")
    plot(x_a(1,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\xi [m]")
subplot(4,1,2)
    hold on; grid on;
    title("\xiDot vs. k")
    plot(x_a(2,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\xiDot [m/s]")
subplot(4,1,3)
    hold on; grid on;
    title("\eta vs. k")
    plot(x_a(3,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\eta [m]")
subplot(4,1,4)
    hold on; grid on;
    title("\etaDot vs. k")
    plot(x_a(4,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\etaDot [m/s]")

% Plot aircraft A's +2sigma bound
figure
sgtitle("$\hat{x}_A(k)$ upper $2\sigma$ bound vs. k", 'Interpreter', 'latex')
subplot(4,1,1)
    hold on; grid on;
    title("2\sigma_\xi vs. k")
    plot(twoSig_a(1,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_\xi [m]")
subplot(4,1,2)
    hold on; grid on;
    title("2\sigma_{\xiDot} vs. k")
    plot(twoSig_a(2,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_{\xiDot} [m/s]")
subplot(4,1,3)
    hold on; grid on;
    title("2\sigma_\eta vs. k")
    plot(twoSig_a(3,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_\eta [m]")
subplot(4,1,4)
    hold on; grid on;
    title("2\sigma_{\etaDot} vs. k")
    plot(twoSig_a(4,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_{\etaDot} [m/s]")

% Plot aircraft B's state
figure
sgtitle("$\hat{x}_B(k)$ vs. k", 'Interpreter', 'latex')
subplot(4,1,1)
    hold on; grid on;
    title("\xi vs. k")
    plot(x_b(1,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\xi [m]")
subplot(4,1,2)
    hold on; grid on;
    title("\xiDot vs. k")
    plot(x_b(2,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\xiDot [m/s]")
subplot(4,1,3)
    hold on; grid on;
    title("\eta vs. k")
    plot(x_b(3,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\eta [m]")
subplot(4,1,4)
    hold on; grid on;
    title("\etaDot vs. k")
    plot(x_b(4,:))
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("\etaDot [m/s]")

% Plot aircraft B's +2sigma bound
figure
sgtitle("$\hat{x}_B(k)$ upper $2\sigma$ bound vs. k", 'Interpreter', 'latex')
subplot(4,1,1)
    hold on; grid on;
    title("2\sigma_\xi vs. k")
    plot(twoSig_b(1,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_\xi [m]")
subplot(4,1,2)
    hold on; grid on;
    title("2\sigma_{\xiDot} vs. k")
    plot(twoSig_b(2,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_{\xiDot} [m/s]")
subplot(4,1,3)
    hold on; grid on;
    title("2\sigma_\eta vs. k")
    plot(twoSig_b(3,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_\eta [m]")
subplot(4,1,4)
    hold on; grid on;
    title("2\sigma_{\etaDot} vs. k")
    plot(twoSig_b(4,:), 'r--')
    xlabel("k (\DeltaT = 0.5 sec)"); ylabel("2\sigma_{\etaDot} [m/s]")

%% Extra
% Plot both aircraft trajectories
figure
hold on; grid on; axis equal
title("Aircraft Propagated Trajectories")
traj_a = plot(x_a(1,:), x_a(3,:), 'b-');
data_a = plot(yaugHist(1,:), yaugHist(2,:), 'c--');
start = plot(x_a(1,1), x_a(3,1), 'g.', 'MarkerSize', 15);
stop = plot(x_a(1,end), x_a(3,end), 'r.', 'MarkerSize', 15);
traj_b = plot(x_b(1,:), x_b(3,:), 'r-');
data_b = plot(x_bData(1,:), x_bData(2,:), '--', 'Color', [1 0.5 0.1]);
plot(x_b(1,1), x_b(3,1), 'g.', 'MarkerSize', 15);
plot(x_b(1,end), x_b(3,end), 'r.', 'MarkerSize', 15);
xlim([-4000 6000]); ylim([-4000 6000])
xlabel("\xi [m]"), ylabel("\eta [m]")
legend([traj_a, data_a, traj_b, data_b, start, stop], ...
       ["Estimated Aircraft A Trajectory", "Measured Aircraft A Trajectory", "Estimated Aircraft B Trajectory", ...
        "Extracted Aircraft B Trajectory", "Start of trajectory", "End of Trajectory"], ...
        'location', 'bestoutside')




