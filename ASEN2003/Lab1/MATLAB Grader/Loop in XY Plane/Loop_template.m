%% User Inputs
clc; clear; close all;

prop = 0.8;
xyz_start = [20 0 20];
vel_in = [50 0 0];
r = 40;
n = 100;
h0 = 125;

%%%%%%%%%%%% MATLAB Grader Inputs - DO NOT CHANGE %%%%%%%%%%%
path_ref = NaN;
theta_vec_ref = NaN;
theta_in_ref = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solution
[theta_in, path, distance, theta_vec, vel_out, g_s] = loop(vel_in, xyz_start, n, prop, r, h0, theta_in_ref, path_ref, theta_vec_ref);

%% Visualization

% Plot path
figure
sgtitle("Paths for r=" + r + ", prop=" + prop)

subplot(3,1,1)
plot(path(:,1), 'b','linewidth', 1.5)
title('Roller Coaster X-Dir Path')
grid on

subplot(3,1,2)
plot(path(:,2), 'b','linewidth', 1.5)
title('Roller Coaster Y-Dir Path')
grid on

subplot(3,1,3)
plot(path(:,3), 'b','linewidth', 1.5)
title('Roller Coaster Z-Dir Path')
grid on

figure
hold on
scatter3(path(1,1), path(1,2), path(1,3));
plot3(path(:,1),path(:,2),path(:,3))
view([30 35])
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')

% Plot gs 
figure
grid on
plot(distance, g_s(:,3), 'b', 'linewidth', 1.5)
title("Vertical Gs for r=" + r + ", prop=" + prop)
ylabel('g''s')
grid on
