%% User Inputs
clc; clear; close all;

xyz_start = [20 0 20];
vel_in = [32 0 32];
n = 100;
h0 = 125;
a = -9.81;
t = (-vel_in(3) - sqrt(vel_in(3)^2 - 2*a*(xyz_start(3)-xyz_start(3))))/(a);
d = vel_in(1)*t;

%%%%%%%%%%%% MATLAB Grader Inputs - DO NOT CHANGE %%%%%%%%%%%
path_ref = NaN;
theta_vec_ref = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solution
[path, distance, vel_out, g_s] = parabola(vel_in, xyz_start, n, a, d, h0, path_ref, theta_vec_ref);

%% Visualization

% Plot path
figure
sgtitle("Paths for a=" + a + ", d=" + d)

subplot(3,1,1)
plot(path(:,1), 'b','linewidth', 1.5)
title('Roller Coaster X-Dir Path - Reference Solution and Your Solution')
grid on

subplot(3,1,2)
plot(path(:,2), 'b','linewidth', 1.5)
title('Roller Coaster Y-Dir Path - Reference Solution and Your Solution')
grid on

subplot(3,1,3)
plot(path(:,3), 'b','linewidth', 1.5)
title('Roller Coaster Z-Dir Path - Reference Solution and Your Solution')
grid on

figure
hold on;
plot3(path(:,1), path(:,2), path(:,3));
view([0,0]);

% Plot g-forces 
figure
grid on
plot(distance, g_s(:,3), 'b', 'linewidth', 1.5)
title("Vertical G''s for a=" + a + ", d=" + d)
ylabel('g''s)')
grid on