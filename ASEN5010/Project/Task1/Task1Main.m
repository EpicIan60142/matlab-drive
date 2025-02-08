%% ASEN 5010 Task 1 main script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath('..\..\Utilities\')

R_Mars = 3396.19; % km

[marsX, marsY, marsZ] = sphere(100);
marsX = R_Mars*marsX;
marsY = R_Mars*marsY;
marsZ = R_Mars*marsZ;

w_LMO = [0; 0; 0.000884797]; % rad/s, O frame coords
EA_LMO = deg2rad([20; 30; 60]); % Omega, i, theta
h_LMO = 400; % km
radius_LMO = R_Mars + h_LMO; % km
x0_LMO = [radius_LMO; EA_LMO; w_LMO];

w_GMO = [0; 0; 0.0000709003]; % rad/s, O frame coords
EA_GMO = deg2rad([0; 0.0000000; 250]); % Omega, i, theta
h_GMO = 17028.01; % km
radius_GMO = R_Mars + h_GMO; % km
x0_GMO = [radius_GMO; EA_GMO; w_GMO];

t0 = 0; % sec
dt = 1; % sec
tf = 6500; % sec

%% Propagate orbits
out_LMO = RK4_Orbit(x0_LMO, t0, dt, tf);
out_GMO = RK4_Orbit(x0_GMO, t0, dt, tf);

%% Extract answers to text files
t_LMO = 450;
t_GMO = 1150;


LMO_ans = out_LMO(out_LMO(:,1) == t_LMO,:);
r_LMO = LMO_ans(2:4);
v_LMO = LMO_ans(5:7);

f1 = fopen("LMO_1.txt", "w");
ans_LMO_1 = fprintf(f1, "%.1f %.1f %.1f", r_LMO(1), r_LMO(2), r_LMO(3));
fclose(f1);

f2 = fopen("LMO_2.txt", "w");
ans_LMO_2 = fprintf(f2, "%.3f %.3f %.3f", v_LMO(1), v_LMO(2), v_LMO(3));
fclose(f2);


GMO_ans = out_GMO(out_GMO(:,1) == t_GMO,:);
r_GMO = GMO_ans(2:4);
v_GMO = GMO_ans(5:7);

f3 = fopen("GMO_1.txt", "w");
ans_GMO_1 = fprintf(f3, "%.1f %.1f %.1f %.1f", r_GMO(1), r_GMO(2), r_GMO(3));
fclose(f3);

f4 = fopen("GMO_2.txt", "w");
ans_GMO_2 = fprintf(f4, "%.3f %.3f %.3f", v_GMO(1), v_GMO(2), v_GMO(3));
fclose(f4);


%% Plot for error checking
fig = figure;

% ax2 = axes();
% I2 = imread("marsStars.jpg");
% imshow(I2, 'parent', ax2)

ax1 = axes();
title("Orbit Simulation", 'Color', 'w')
hold on
grid on
axis equal
plot3(out_LMO(:,2), out_LMO(:,3), out_LMO(:,4), 'LineWidth', 3)
plot3(out_GMO(:,2), out_GMO(:,3), out_GMO(:,4), 'LineWidth', 3)

mars = surf(marsX, marsY, -marsZ); % Need to flip the sphere for image to map properly
I = imread("marsSurface.jpg");
set(mars,'FaceColor','texturemap','cdata',I,'edgecolor','none');

% set(gca, 'Color', 'none', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', 'GridColor', 'w')

view([30 35])

