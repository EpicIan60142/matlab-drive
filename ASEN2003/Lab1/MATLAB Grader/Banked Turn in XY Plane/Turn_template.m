clc; clear; close all;

%% User Inputs
r = 50;
bank_angle = 50*pi/180;
turn_dir = -1;
pos_in = [10 0 0];
vel_in = [49.5227 0 0];
h0 = 125;
n = 100;

%%%%%%%%%%%% MATLAB Grader Inputs - DO NOT CHANGE %%%%%%%%%%%
path_ref = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solution
[path, distance, vel_out, g_forces] = banked_turn(vel_in, pos_in, n, r, turn_dir, bank_angle, h0, path_ref);

%% Visualization

% Plot path
figure
sgtitle("Paths for r=" + r + ", turndir=" + turn_dir)

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
view([30, 35])

% Plot gs 
figure
grid on
hold on
plot(distance, g_forces(:,3), 'linewidth', 1.5)
plot(distance, g_forces(:,2), 'linewidth', 1.5)
plot(distance, g_forces(:,1), 'linewidth', 1.5)
title("Gs vs Distance along Element for r=" + r + ", turndir=" + turn_dir)
legend('Vertical Gs', 'Lateral Gs', 'Front/Back Gs')
ylabel('g''s)')
xlabel('Distance(m)')
grid on