%% ASEN 5014 Project 1 Main Script
% By: Ian Faber
% Partner: Gabriel Agostine

%% Housekeeping
clc; clear; close all;

%% Setup
format shortE

digits = 16; % Number of digits to keep for numerical precision

muEarth = 398600; % km^3/s^2
rOrbit = 6778; % km - 400 km LEO orbit

syms n lambda

A = [
        0       0       0       1       0   0;
        0       0       0       0       1   0;
        0       0       0       0       0   1;
        3*n^2   0       0       0       2*n 0;
        0       0       0       -2*n    0   0;
        0       0       -n^2    0       0   0
    ];

B = [
        0 0 0;
        0 0 0;
        0 0 0;
        1 0 0;
        0 1 0;
        0 0 1
    ];

C = [
        1 0 0 0 0 0;
        0 1 0 0 0 0;
        0 0 1 0 0 0
    ];

D = [
        0 0 0;
        0 0 0;
        0 0 0
    ];

%% Part 1
fprintf("--- Part 1 ---\n")
% x = [x; y; z; xDot; yDot; zDot], y = [x; y; z], u = [ux; uy; uz]
charEq = charpoly(A, lambda) % Get characteristic equation

n = sqrt(muEarth/(rOrbit^3));

    % Implement plant
A = [
        0       0       0       1       0   0;
        0       0       0       0       1   0;
        0       0       0       0       0   1;
        3*n^2   0       0       0       2*n 0;
        0       0       0       -2*n    0   0;
        0       0       -n^2    0       0   0
    ];

[G, Lambda] = eig(A); % Eigenvectors and eigenvalues

    % Fix numerical precision
G = round(G,digits);
Lambda = round(Lambda, digits)

if det(G) == 0 % G is singular, need more eigenvectors
    [V, J] = jordan(A); % Extended eigenvectors and jordan form
    X1 = V(:,1);
    X2 = V(:,2);
    U3 = real(V(:,3));
    V3 = imag(V(:,3));
    U4 = real(V(:,5));
    V4 = imag(V(:,5));
else
    X1 = G(:,1);
    X2 = G(:,2);
    U3 = real(G(:,3));
    V3 = imag(G(:,3));
    U4 = real(G(:,5));
    V4 = imag(G(:,5));
end

    % Eigenspaces and real modal spaces
eigenSpaces = V
modeVecs = [X1, X2, U3, V3, U4, V4]

realModalForm = (modeVecs^-1)*A*modeVecs;

realModalForm = round(realModalForm, digits)

    % Simulate modes
dt = 1; 
t_modes = 0:dt:10000;
phi = @(t) expm(A*t);
for k = 1:4
    switch k
        case 1
            x0 = X1;
        case 2
            x0 = X2;
        case 3
            x0 = U3;
        case 4
            x0 = U4;
        otherwise
            x0 = zeros(6,1); % Should never reach this
    end
    % x0 = modeVecs(:,k);
    x = zeros(6, length(t_modes));
    for kk = 1:length(t_modes)
        x(:,kk) = phi(t_modes(kk))*x0;
    end

    titleText = sprintf("System open loop response along mode %.0f_{ol} vs. time:\nx_0 = [%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]^T", k, x0);
    figure
    sgtitle(titleText, 'FontSize', 11)
    subplot(3,2,1)
        title("x vs. time")
        hold on; grid on;
        plot(t_modes, x(1,:));
        if k == 2
            ylim([3.5 4.5]);
        end
        xlabel("Time [sec]"); ylabel("x [km]")
    subplot(3,2,2)
        title("xDot vs. time")
        hold on; grid on;
        plot(t_modes, x(4,:));
        if k == 2
            ylim([-0.5 0.5]);
        end
        xlabel("Time [sec]"); ylabel("xDot [km/s]")
    subplot(3,2,3)
        title("y vs. time")
        hold on; grid on;
        plot(t_modes, x(2,:));
        if k == 1
            ylim([-7.5e-3 -6.5e-3])
        end
        xlabel("Time [sec]"); ylabel("y [km]")
    subplot(3,2,4)
        title("yDot vs. time")
        hold on; grid on;
        plot(t_modes, x(5,:));
        if k == 2
            ylim([-7.5e-3 -6.5e-3])
        end
        xlabel("Time [sec]"); ylabel("yDot [km/s]")
    subplot(3,2,5)
        title("z vs. time")
        hold on; grid on;
        plot(t_modes, x(3,:));
        xlabel("Time [sec]"); ylabel("z [km]")
    subplot(3,2,6)
        title("zDot vs. time")
        hold on; grid on;
        plot(t_modes, x(6,:));
        xlabel("Time [sec]"); ylabel("zDot [km/s]")
    
    titleText = sprintf("System open loop response along mode %.0f in x-y-z space:\nx_0 = [%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]^T", k, x0);
    figure
    hold on; grid on;
    title(titleText)

    trajectory = plot3(x(1,:), x(2,:), x(3,:), 'b.', 'MarkerSize', 5);
    start = plot3(x(1,1), x(2,1), x(3,1), 'g.', 'MarkerSize', 15);
    stop = plot3(x(1,end), x(2,end), x(3,end), 'r.', 'MarkerSize', 15);

    xlabel("x [km]"); ylabel("y [km]"); zlabel("z [km]")
    legend([start, trajectory, stop], ["Initial Condition", "Trajectory", "End of Trajectory"], "Location", "best")

    switch k
        case 1
            xlim([-0.5 0.5]); ylim([-7.5e-3 -6.3e-3]); zlim([-0.5 0.5]);
            view([30 35]);
        case 2
            xlim([3.5 4.5]);
            view([30 35]);
        case 3
            xlim([-2 2]); ylim([-3.5 3.5])
            view([0 90]);
        case 4
            zlim([-0.75 0.75]);
            view([30 35]);
    end
end

% Since zero eigenvalues, plant is Lyapunov stable in all states except y,
% which has a secular drift

%% Part 2
fprintf("--- Part 2 ---\n")
P = [B, A*B, A^2*B, A^3*B, A^4*B, A^5*B];
reachRank = rank(P)

orthBasis_reachable = eye(size(A))

% orthBasis = orth(P)

% [orthBasis, ~] = qr(P)

T = 2*pi/n;
% t01 = T/10; % t1 - t0
t01 = T;

    % Perturb A by a tiny amount for lyap
Aperturb = A - 1e-6*eye(size(A));%1e-10*diag([1, 1, 1, 0, 0, 0]);

    % Build Q
Q = expm(-Aperturb*t01)*(B*B')*expm(-Aperturb'*t01) - B*B';
G = lyap(Aperturb, Q)

    % Find minimum energy control energy and v vectors
energy_min = [];
v = [];
for k = 1:size(orthBasis_reachable, 2)
    zeta = zeros(size(orthBasis_reachable(:,k))) - orthBasis_reachable(:,k);
    v(:,k) = (G^-1)*zeta;
    energy_min = [energy_min; zeta'*(G^-1)*zeta];
end
energy_min

    % Sanity check
sys_ol = ss(A,B,eye(6),zeros(size(B)));

energy_ol = [];
t_ol = 0:dt:t01;
for k = 1:size(v,2)
    x0 = orthBasis_reachable(:,k);
    u_ol = [];
    for kk = 1:length(t_ol)
        u_ol = [u_ol, B'*expm(-A'*t_ol(kk))*v(:,k)];
    end
    
    x = lsim(sys_ol, u_ol, t_ol, x0);
    x = x';

    titleText = sprintf("System minimum energy control perturbation response vs. time:\n\\deltax_0 = [%.0f, %.0f, %.0f, %.0f, %.0f, %.0f]^T", x0);
    figure
    sgtitle(titleText, 'FontSize', 11)
    subplot(3,3,1)
        title("x vs. time")
        hold on; grid on;
        plot(t_ol, x(1,:));
        xlabel("Time [sec]"); ylabel("x [km]")
    subplot(3,3,2)
        title("xDot vs. time")
        hold on; grid on;
        plot(t_ol, x(4,:));
        xlabel("Time [sec]"); ylabel("xDot [km/s]")
    subplot(3,3,3)
        title("u_1 vs. time")
        hold on; grid on;
        plot(t_ol, u_ol(1,:));
        xlabel("Time [sec]"); ylabel("u_1 [km/s^2]")
    subplot(3,3,4)
        title("y vs. time")
        hold on; grid on;
        plot(t_ol, x(2,:));
        xlabel("Time [sec]"); ylabel("y [km]")
    subplot(3,3,5)
        title("yDot vs. time")
        hold on; grid on;
        plot(t_ol, x(5,:));
        xlabel("Time [sec]"); ylabel("yDot [km/s]")
    subplot(3,3,6)
        title("u_2 vs. time")
        hold on; grid on;
        plot(t_ol, u_ol(2,:));
        xlabel("Time [sec]"); ylabel("u_2 [km/s^2]")
    subplot(3,3,7)
        title("z vs. time")
        hold on; grid on;
        plot(t_ol, x(3,:));
        xlabel("Time [sec]"); ylabel("z [km]")
    subplot(3,3,8)
        title("zDot vs. time")
        hold on; grid on;
        plot(t_ol, x(6,:));
        xlabel("Time [sec]"); ylabel("zDot [km/s]")
    subplot(3,3,9)
        title("u_3 vs. time")
        hold on; grid on;
        plot(t_ol, u_ol(3,:));
        xlabel("Time [sec]"); ylabel("u_3 [km/s^2]")
end



%% Part 3
fprintf("--- Part 3 ---\n")

    % Choose time constant
steadytime = linspace(600, 1500, 6);
tau = steadytime/5;
% tau = T/10; % On the timescale of an orbit

    % Create system poles
eigVals_cl = -1./tau
% eigVals_cl = [-1/tau, -1/tau, -2/tau, -2/tau, -3/tau, -3/tau]


%% Part 4
fprintf("--- Part 4 ---\n")

    % Design feedback controller
K = place(A, B, eigVals_cl)

Acl = A-(B*K);

    % Find new eigenvectors and closed loop modal spaces
[v_cl, eig_cl] = eig(Acl)

phi_cl = @(t)expm(Acl*t);

    % Simulate closed loop system and find closed loop energy
t_cl = 0:dt:t01;
t_clol = t_cl(t_cl < Inf);
energy_cl = [];
for k = 1:6
        % Simulate mode perturbation
    x0 = orthBasis_reachable(:,k);
    x = zeros(6, length(t_cl));
    u_cl = zeros(3, length(t_cl));
    u_ol = zeros(3, length(t_clol));
    eng = [];
    for kk = 1:length(t_cl)
        x(:,kk) = phi_cl(t_cl(kk))*x0;
        u_cl(:,kk) = -K*x(:,kk); % Closed loop control signal
        eng = [eng; u_cl(:,kk)'*u_cl(:,kk)];  
    end
    for kk = 1:length(t_clol)
        u_ol(:,kk) = B'*expm(-A'*t_clol(kk))*v(:,k); % Open loop control signal
    end
    energy_cl = [energy_cl; sum(eng)];

        % Make time domain plots
    titleText = sprintf("System closed loop perturbation response vs. time:\n\\deltax_0 = [%.0f, %.0f, %.0f, %.0f, %.0f, %.0f]^T", x0);
    figure
    sgtitle(titleText, 'FontSize', 11)
    subplot(3,3,1)
        title("x vs. time")
        hold on; grid on;
        plot(t_cl, x(1,:));
        xlabel("Time [sec]"); ylabel("x [km]")
    subplot(3,3,2)
        title("xDot vs. time")
        hold on; grid on;
        plot(t_cl, x(4,:));
        xlabel("Time [sec]"); ylabel("xDot [km/s]")
    subplot(3,3,3)
        title("u_1 vs. time")
        hold on; grid on;
        plot(t_cl, u_cl(1,:));
        xlabel("Time [sec]"); ylabel("u_1 [km/s^2]")
    subplot(3,3,4)
        title("y vs. time")
        hold on; grid on;
        plot(t_cl, x(2,:));
        xlabel("Time [sec]"); ylabel("y [km]")
    subplot(3,3,5)
        title("yDot vs. time")
        hold on; grid on;
        plot(t_cl, x(5,:));
        xlabel("Time [sec]"); ylabel("yDot [km/s]")
    subplot(3,3,6)
        title("u_2 vs. time")
        hold on; grid on;
        plot(t_cl, u_cl(2,:));
        xlabel("Time [sec]"); ylabel("u_2 [km/s^2]")
    subplot(3,3,7)
        title("z vs. time")
        hold on; grid on;
        plot(t_cl, x(3,:));
        xlabel("Time [sec]"); ylabel("z [km]")
    subplot(3,3,8)
        title("zDot vs. time")
        hold on; grid on;
        plot(t_cl, x(6,:));
        xlabel("Time [sec]"); ylabel("zDot [km/s]")
    subplot(3,3,9)
        title("u_3 vs. time")
        hold on; grid on;
        plot(t_cl, u_cl(3,:));
        xlabel("Time [sec]"); ylabel("u_3 [km/s^2]")
    
        % Make x-y-z plots
    titleText = sprintf("System closed loop perturbation response in x-y-z space:\n\\deltax_0 = [%.0f, %.0f, %.0f, %.0f, %.0f, %.0f]^T", x0);
    figure
    hold on; grid on;
    title(titleText)
    trajectory = plot3(x(1,:), x(2,:), x(3,:), 'b.', 'MarkerSize', 5);
    start = plot3(x(1,1), x(2,1), x(3,1), 'g.', 'MarkerSize', 15);
    stop = plot3(x(1,end), x(2,end), x(3,end), 'r.', 'MarkerSize', 15);
    xlabel("x [km]"); ylabel("y [km]"); zlabel("z [km]")
    legend([start, trajectory, stop], ["Initial Condition", "Trajectory", "End of Trajectory"], "Location", "best")
    view([30 35])

        % Make control signal plots
    titleText = sprintf("Control signal comparison:\n\\deltax_0 = [%.0f, %.0f, %.0f, %.0f, %.0f, %.0f]^T", x0);
    figure; t = tiledlayout(3,1);
    title(t, titleText, 'FontSize', 11); 
    nexttile; 
        hold on; grid on;
        title("u_1 vs. time")
        minPlot = plot(t_clol, u_ol(1,:),'b');
        clPlot = plot(t_cl, u_cl(1,:),'r--');
        xlabel("Time [sec]"); ylabel("u_1 [km/s^2]")
    nexttile; 
        hold on; grid on;
        title("u_2 vs. time")
        plot(t_clol, u_ol(2,:),'b')
        plot(t_cl, u_cl(2,:),'r--');
        xlabel("Time [sec]"); ylabel("u_2 [km/s^2]")
    nexttile; 
        hold on; grid on;
        title("u_3 vs. time")
        plot(t_clol, u_ol(3,:),'b')
        plot(t_cl, u_cl(3,:),'r--');
        xlabel("Time [sec]"); ylabel("u_3 [km/s^2]")
    legend([minPlot, clPlot], ["Minimum Energy Control", "State Feedback Control"], 'location', 'bestoutside')

end
energy_cl

    % Compare closed loop and open loop energies
energy_diff = energy_cl - energy_min

%% Part 5
fprintf("--- Part 5 ---\n");

    % Compute input gain matrix
F = (C*((-A+(B*K))^-1)*B)^-1

    % xDot = (A-BK)x + BFr
    % y = (C-DK)x + DFr, D = zeros so y = Cx as before
B_F = B*F;
Cplot = eye(size(A));
Dplot = zeros(size(B));

    % Plot individual step responses
sys_cl = ss(Acl, B, Cplot, Dplot, 'StateName', {'x','y','z','vx','vy','vz'}, 'InputName', {'u1','u2','u3'}, 'OutputName', {'x [km]','y [km]','z [km]','vx [km/s]','vy [km/s]','vz [km/s]'});
sys_F = ss(Acl, B_F, Cplot, Dplot, 'StateName', {'x','y','z','vx','vy','vz'}, 'InputName', {'u1','u2','u3'}, 'OutputName', {'x [km]','y [km]','z [km]','vx [km/s]','vy [km/s]','vz [km/s]'});

stepOpt = timeoptions("cstprefs");
stepOpt.YLabel.String = "State Response";
stepOpt.Grid = 'on';
stepOpt.OutputLabels.FontSize = 7.5;

stepOpt_F = timeoptions("cstprefs");
stepOpt_F.YLim = {[-1 1], [-1 1], [-1 1], [-3e-3 3e-3], [-3e-3 3e-3], [-3e-3 3e-3]};
stepOpt_F.YLabel.String = "State Response";
stepOpt_F.Grid = 'on';
stepOpt_F.OutputLabels.FontSize = 7.5;

figure
stepplot(sys_cl, t01, stepOpt)
title("System closed loop step response to unit inputs - no F matrix")

figure
stepplot(sys_F, t01, stepOpt_F)
title("System closed loop step response to unit inputs - F matrix included")

%% Reset format to default
format default





