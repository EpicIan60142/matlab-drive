%% ASEN 6014 HW 3 Problem 4 Main Script
% By: Ian Faber, 11/16/2025

%% Housekeeping
clc; clear; close all;

%% Setup
    % Planetary constants
pConst.mu = 398600.4415; % km^3/s^2
pConst.Ri = 6378; % km
pConst.J2 = 1.08264e-3;

    % Chief initial orbit elements
oe_c0.mu = pConst.mu;
oe_c0.a = 7555; % km
oe_c0.e = 0.05;
oe_c0.i = deg2rad(48); % deg -> rad
oe_c0.RAAN = deg2rad(20); % deg -> rad
oe_c0.argPeri = deg2rad(10); % deg -> rad
M_c = deg2rad(20); % deg -> rad
f_c = convE2f(convM2E(M_c, oe_c0.e, false), oe_c0.e);
oe_c0.f = f_c;

    % Desired deputy orbit element differences - from slide 57
doe_r.da = -0.00192995; % km
doe_r.de = 0.000576727; 
doe_r.di = deg2rad(0.006);
doe_r.dRAAN = deg2rad(0);
doe_r.dargPeri = deg2rad(0);
doe_r.dM0 = deg2rad(0.1);

    % Initial deputy minus desired orbit element differences
doe.da = 0.1; % km
doe.de = 0;
doe.di = deg2rad(-0.05);
doe.dRAAN = deg2rad(0.01);
doe.dargPeri = deg2rad(0);
doe.dM0 = deg2rad(0);

    % Initial deputy orbit elements
oe_d = oe_c0;
oe_d.a = oe_c0.a + doe.da + doe_r.da;
oe_d.e = oe_c0.e + doe.de + doe_r.de;
oe_d.i = oe_c0.i + doe.di + doe_r.di;
oe_d.RAAN = oe_c0.RAAN + doe.dRAAN + doe_r.dRAAN;
oe_d.argPeri = oe_c0.argPeri + doe.dargPeri + doe_r.dargPeri;
M_d = M_c + doe.dM0 + doe_r.dM0;
f_d = convE2f(convM2E(M_d, oe_d.e, false), oe_d.e);
oe_d.f = f_d;

    % ode45 settings
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);

%% Run controller
    % Setup initial state vector
cart_c = convClassicOE2Cart(oe_c0);
X0_c = [cart_c.rVec; cart_c.vVec];

cart_d = convClassicOE2Cart(oe_d);
X0_d = [cart_d.rVec; cart_d.vVec];

X0 = [X0_c; X0_d];

    % Define tspan for multiple orbits
nOrbits = 4;
T = 2*pi*sqrt((oe_c0.a)^3/pConst.mu);
tspan = 0:10:nOrbits*T;

    % Run controller
[t, X, u, doe, oe_d] = impulsiveFeedbackControl(X0, doe_r, tspan, pConst, opt);
% [t, X] = ode15s(@(t,X)impulsiveFeedbackControl(t,X,doe_r,M_d,pConst), tspan, X0, opt);
% [~, u, doe, oe_d, oe_c] = cellfun(@(t,X)impulsiveFeedbackControl(t,X.',doe_r,M_d,pConst), num2cell(t), num2cell(X,2), 'uni', 0);
% u = cellfun(@(x)x',u,'uni',0);
% doe = cellfun(@(x)x',doe,'uni',0);
tMan = cell2mat(oe_d(:,2));
tMan = [tMan; t(end)];
oe_d = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_d(:,1), 'uni', 0);
% oe_c = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_c(:,1), 'uni', 0);


%% Plot results
    % Extract states and control
X_c = X(:,1:6);
X_d = X(:,7:12);
% u = cell2mat(u);
% doe = cell2mat(doe);
doe = [doe; doe(end,:)];
oe_d = cell2mat(oe_d);
oe_d = [oe_d; oe_d(end,:)];
% oe_c = cell2mat(oe_c);
% oe_c = [oe_c; oe_c(end,:)];

    % Convert to Hill Frame
X_d_Hill = zeros(size(X_d));
for k = 1:length(t)
    X_d_Hill(k,:) = convDeputyN2H(X_c(k,:)', X_d(k,:)', pConst)';
end

    % Plot
markerSize = 30;
figure; axis equal
hold on; grid on;
title("Impulsive Relative Motion Control")
scatter3(X_d_Hill(1,1), X_d_Hill(1,2), X_d_Hill(1,3), markerSize, 'g', 'filled');
plot3(X_d_Hill(:,1), X_d_Hill(:,2), X_d_Hill(:,3), 'b-')
scatter3(X_d_Hill(end,1), X_d_Hill(end,2), X_d_Hill(end,3), markerSize, 'r', 'filled')
% scatter3(0,0,0,markerSize,'k','*')
xlabel("Radial [km]"); ylabel("Along-Track [km]"); zlabel("Cross-Track [km]")
view([180+30,35]);
legend("Start", "Trajectory", "End")%, "Chief")

figure; tl = tiledlayout(3,1); ax = [];
title(tl, "Control Effort vs. Time")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    stairs(t,u(:,1),'b-');
    xlabel("Time [sec]"); ylabel("\Delta v_r [km/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    stairs(t,u(:,2),'b-');
    xlabel("Time [sec]"); ylabel("\Delta v_\theta [km/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    stairs(t,u(:,3),'b-');
    xlabel("Time [sec]"); ylabel("\Delta v_h [km/s]");
linkaxes(ax, 'x');

figure; tl = tiledlayout(3,1); %ax = [];
title(tl, "Velocity vs. Time")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t,X_d_Hill(:,4),'b-');
    xlabel("Time [sec]"); ylabel("v_x [km/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t,X_d_Hill(:,5),'b-');
    xlabel("Time [sec]"); ylabel("v_y [km/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t,X_d_Hill(:,6),'b-');
    xlabel("Time [sec]"); ylabel("v_z [km/s]");
linkaxes(ax, 'x');

figure; tl = tiledlayout(3,2);
title(tl, "Deputy Elements vs. Time")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    stairs(tMan,oe_d(:,1), 'b-');
    xlabel("Time [sec]"); ylabel("a_d [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    stairs(tMan,oe_d(:,2),'b-');
    xlabel("Time [sec]"); ylabel("e_d")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    stairs(tMan,rad2deg(oe_d(:,3)),'b-');
    xlabel("Time [sec]"); ylabel("i_d [deg]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    stairs(tMan,rad2deg(oe_d(:,4)),'b-');
    xlabel("Time [sec]"); ylabel("\Omega_d [deg]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    stairs(tMan,rad2deg(oe_d(:,5)),'b-');
    xlabel("Time [sec]"); ylabel("\omega_d [deg]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(tMan,rad2deg(oe_d(:,6)),'b-');
    xlabel("Time [sec]"); ylabel("M_d [deg]")
linkaxes(ax, 'x');

figure; tl = tiledlayout('flow'); %ax = [];
title(tl, "Difference in Desired and Actual Deputy Elements vs. Time")
try doe_r.da;
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        stairs(tMan,doe(:,1),'b-');
        plot(t,zeros(size(t)), 'r--')
        xlabel("Time [sec]"); ylabel("\Delta a [km]")
end
try doe_r.de;
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        stairs(tMan,doe(:,2),'b-');
        plot(t,zeros(size(t)), 'r--')
        xlabel("Time [sec]"); ylabel("\Delta e")
end
try doe_r.di;
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        stairs(tMan,rad2deg(doe(:,3)),'b-');
        plot(t,zeros(size(t)), 'r--')
        xlabel("Time [sec]"); ylabel("\Delta i [deg]")
end
try doe_r.dRAAN;
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        stairs(tMan,rad2deg(doe(:,4)),'b-');
        plot(t,zeros(size(t)), 'r--')
        xlabel("Time [sec]"); ylabel("\Delta \Omega [deg]")
end
try doe_r.dargPeri;
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        stairs(tMan,rad2deg(doe(:,5)),'b-');
        plot(t,zeros(size(t)), 'r--')
        xlabel("Time [sec]"); ylabel("\Delta \omega [deg]")
end
try doe_r.dM0;
    nt = nexttile; ax = [ax; nt];
        hold on; grid on;
        plot(tMan,rad2deg(doe(:,6)),'b-');
        plot(t,zeros(size(t)), 'r--')
        xlabel("Time [sec]"); ylabel("\Delta M [deg]")
end
linkaxes(ax, 'x');

% return;
    % Animate
frames = 10;
figure;
for k = 1:frames:length(t)
    clf;
    hold on; grid on; axis equal
    title(sprintf("Deputy Hill Trajectory at t = %.3f sec", t(k)));
    cutoff = find(t >= 2*T, 1, 'first');
    if k <= cutoff
        start = 1;
        plot3(X_d_Hill(1,1), X_d_Hill(1,2), X_d_Hill(1,3), 'g.', 'MarkerSize', markerSize);
    else
        start = floor(k - cutoff);
    end
    plot3(X_d_Hill(start:k,1), X_d_Hill(start:k,2), X_d_Hill(start:k,3), 'b-');
    if k <= length(t) - frames
        plot3(X_d_Hill(k,1), X_d_Hill(k,2), X_d_Hill(k,3), 'k.', 'MarkerSize', markerSize);
    end
    plot3(0,0,0,'k*','MarkerSize',10);
    if k >= length(t) - frames
       plot3(X_d_Hill(end,1), X_d_Hill(end,2), X_d_Hill(end,3), 'r.', 'MarkerSize', markerSize); 
    end
    xlabel("Radial [km]"); ylabel("Along-Track [km]"); zlabel("Cross-Track [km]");
    view([30+k/frames,35]);
    drawnow;
end

