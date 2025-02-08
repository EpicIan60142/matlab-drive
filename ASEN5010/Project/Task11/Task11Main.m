%% ASEN 5010 Task 11 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all 

%% Setup
% Include all of our lovely task functions
addpath('..\..\Utilities')
addpath('..\Task1')
addpath('..\Task2')
addpath('..\Task3')
addpath('..\Task4')
addpath('..\Task5')
addpath('..\Task6')

% Make Mars
R_Mars = 3396.19; % km
[marsX, marsY, marsZ] = sphere(100);
marsX = R_Mars*marsX;
marsY = R_Mars*marsY;
marsZ = R_Mars*marsZ;
marsAngle = 0; % rad

% Initial orbit parameters
h_LMO = 400; % km
h_GMO = 17028.01; % km
radius_LMO = R_Mars + h_LMO; % km
radius_GMO = R_Mars + h_GMO; % km

w_0_LMO = [0; 0; 0.000884797]; % rad/s, O frame coords
EA_0_LMO = deg2rad([20; 30; 60]); % Omega, i, theta

w_0_GMO = [0; 0; 0.0000709003]; % rad/s, O frame coords
EA_0_GMO = deg2rad([0; 0; 250]); % Omega, i, theta

x_0_LMO = [radius_LMO; EA_0_LMO; w_0_LMO];
x_0_GMO = [radius_GMO; EA_0_GMO; w_0_GMO];

% Initial attitude parameters
sigBN_0 = [0.3; -0.4; 0.5]; % unitless
omegBN_0 = deg2rad([1.00; 1.75; -2.20]); % Converted to rad/s from deg/s
u_0 = zeros(3,1); % Nm
I = diag([10, 5, 7.5]); % kgm^2, body coords

x_0_att = {I; sigBN_0; omegBN_0; u_0};

% Controller parameters
K = 1*(1/180); % 0.0056, from zeta requirement
P = 1*(1/6); % 0.1667, from time decay requirement
controlParams = [K; P];

% RK4 params
t0 = 0;
dt = 1;
tf = 6500;
% tf = 90000;

%% Propagate orbits and attitude with full pointing control
out_LMO = RK4_Orbit(x_0_LMO, t0, dt, tf);
out_GMO = RK4_Orbit(x_0_GMO, t0, dt, tf);

[out_Attitude, states_attitude] = RK4_Attitude(x_0_att, t0, dt, tf, out_LMO, out_GMO, controlParams);

%% Extract and process simulation outputs
t = out_Attitude(:,1);
sig = out_Attitude(:,2:4);
w = out_Attitude(:,5:7);
u = out_Attitude(:,8:10);
sigRef = out_Attitude(:,11:13);
wRef = out_Attitude(:,14:16);
states = states_attitude;

idxSun = states == "sun pointing";
idxNadir = states == "nadir pointing";
idxGMO = states == "GMO pointing";

tInt = [300; 2100; 3400; 4400; 5600];

%% Attitude and control plots
figAttitude = figure;

sgtitle("Nano-satellite Attitude Evolution Over Time")

subplot(3,3,1)
hold on
grid on
title("\sigma_1 vs. time")
plot(t, sig(:,1))
plot(t, sigRef(:,1), '--')
xlabel("Time [sec]")
ylabel("\sigma_1")
legend("Current", "Reference", 'location', 'best')

subplot(3,3,2)
hold on
grid on
title("\omega_1 vs. time")
plot(t, w(:,1))
plot(t, wRef(:,1), '--')
xlabel("Time [sec]")
ylabel("\omega_1 [rad/s]")
legend("Current", "Reference", 'location', 'best')

subplot(3,3,3)
hold on
grid on
title("u_1 vs. time")
plot(t, u(:,1))
xlabel("Time [sec]")
ylabel("u_1 [Nm]")

subplot(3,3,4)
hold on
grid on
title("\sigma_2 vs. time")
plot(t, sig(:,2))
plot(t, sigRef(:,2), '--')
xlabel("Time [sec]")
ylabel("\sigma_2")
legend("Current", "Reference", 'location', 'best')

subplot(3,3,5)
hold on
grid on
title("\omega_2 vs. time")
plot(t, w(:,2))
plot(t, wRef(:,2), '--')
xlabel("Time [sec]")
ylabel("\omega_2 [rad/s]")
legend("Current", "Reference", 'location', 'best')

subplot(3,3,6)
hold on
grid on
title("u_2 vs. time")
plot(t, u(:,2))
xlabel("Time [sec]")
ylabel("u_2 [Nm]")

subplot(3,3,7)
hold on
grid on
title("\sigma_3 vs. time")
plot(t, sig(:,3))
plot(t, sigRef(:,3), '--')
xlabel("Time [sec]")
ylabel("\sigma_3")
legend("Current", "Reference", 'location', 'best')

subplot(3,3,8)
hold on
grid on
title("\omega_3 vs. time")
plot(t, w(:,3))
plot(t, wRef(:,3), '--')
xlabel("Time [sec]")
ylabel("\omega_3 [rad/s]")
legend("Current", "Reference", 'location', 'best')

subplot(3,3,9)
hold on
grid on
title("u_3 vs. time")
plot(t, u(:,3))
xlabel("Time [sec]")
ylabel("u_3 [Nm]")

%% Orbit animation
figOrbit = figure('Position',[0 0 1920 1080]);

movieVector = [];
frames = [];
dTime = 50;

for k = 1:dTime:length(t)
    clf
    hold on
    grid on
    axis equal

    % subplot(1,2,1)
    % ax2 = axes();
    % marsStars = imread("marsStars.jpg");
    % imshow(marsStars, 'parent', ax2)
    
    % ax1 = axes();
    titleText = sprintf("ASEN 5010 Capstone: t = %.0f sec, mode = %s", t(k), states(k));
    title(titleText, 'Color', 'k')
    % axis equal

    % Orbits
    nanoOrbit = plot3(out_LMO(:,2), out_LMO(:,3), out_LMO(:,4), 'LineWidth', 2);
    motherOrbit = plot3(out_GMO(:,2), out_GMO(:,3), out_GMO(:,4), 'LineWidth', 2);
    
    % Spacecraft
        % Nano-sat
    nanoNH = calcHN(t(k), out_LMO)';
    nanoDCM = MRP2DCM(out_Attitude(k,2:4))';
    nanoSat = scatter3(out_LMO(k,2), out_LMO(k,3), out_LMO(k,4), 25, 'black', 'filled');
    quiver3(0, 0, 0, nanoNH(1,1), nanoNH(2,1), nanoNH(3,1), 6250, 'k') % o_1 axis
    quiver3(0, 0, 0, nanoNH(1,2), nanoNH(2,2), nanoNH(3,2), 6250, 'k') % o_2 axis
    quiver3(0, 0, 0, nanoNH(1,3), nanoNH(2,3), nanoNH(3,3), 6250, 'k') % o_3 axis
    b_1 = quiver3(out_LMO(k,2), out_LMO(k,3), out_LMO(k,4), nanoDCM(1,1), nanoDCM(2,1), nanoDCM(3,1), 2500, 'r'); % b_1 axis
    negB_1 = quiver3(out_LMO(k,2), out_LMO(k,3), out_LMO(k,4), -nanoDCM(1,1), -nanoDCM(2,1), -nanoDCM(3,1), 2500, 'r:'); % -b_1 axis
    b_2 = quiver3(out_LMO(k,2), out_LMO(k,3), out_LMO(k,4), nanoDCM(1,2), nanoDCM(2,2), nanoDCM(3,2), 2500, 'g'); % b_2 axis
    b_3 = quiver3(out_LMO(k,2), out_LMO(k,3), out_LMO(k,4), nanoDCM(1,3), nanoDCM(2,3), nanoDCM(3,3), 2500, 'b'); % b_3 axis
        % Mothercraft
    GMOSat = scatter3(out_GMO(k,2), out_GMO(k,3), out_GMO(k,4), 25, 'magenta', 'filled');
    GMONH = calcHN(t(k), out_GMO)';
    quiver3(0, 0, 0, GMONH(1,1), GMONH(2,1), GMONH(3,1), 25000, 'k') % o_1 axis
    quiver3(0, 0, 0, GMONH(1,2), GMONH(2,2), GMONH(3,2), 25000, 'k') % o_2 axis
    quiver3(0, 0, 0, GMONH(1,3), GMONH(2,3), GMONH(3,3), 25000, 'k') % o_3 axis

    % Other frames
        % Communication frame
    % RcN = calcRcN(t(k), out_GMO, out_LMO)';
    % quiver3(out_LMO(k,2), out_LMO(k,3), out_LMO(k,4), RcN(1,1), RcN(2,1), RcN(3,1), 10000, 'r') % r_1 axis
    % quiver3(out_LMO(k,2), out_LMO(k,3), out_LMO(k,4), -RcN(1,1), -RcN(2,1), -RcN(3,1), 10000, 'r:') % -r_1 axis
    % quiver3(out_LMO(k,2), out_LMO(k,3), out_LMO(k,4), RcN(1,2), RcN(2,2), RcN(3,2), 10000, 'g') % r_2 axis
    % quiver3(out_LMO(k,2), out_LMO(k,3), out_LMO(k,4), RcN(1,3), RcN(2,3), RcN(3,3), 10000, 'b') % r_3 axis

    % Mars
    marsAngle = marsAngle + dTime*(w_0_GMO(3)); % How much has Mars rotated over the timestep?
    marsXrot = marsX*cos(marsAngle) - marsY*sin(marsAngle);
    marsYrot = marsX*sin(marsAngle) + marsY*cos(marsAngle);
    marsZrot = -marsZ; % Need to flip the sphere for image to map properly
    mars = surf(marsXrot, marsYrot, marsZrot); 
    marsSurface = imread("marsSurface.jpg");
    set(mars,'FaceColor','texturemap','cdata',marsSurface,'edgecolor','none');
    
    % set(gca, 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', 'GridColor', 'k')
    xlim([-21000 21000])
    ylim([-21000 21000])
    zlim([-15000 15000])
    xlabel("$\hat{n}_1$", "Interpreter","latex")
    ylabel("$\hat{n}_2$", "Interpreter","latex")
    zlabel("$\hat{n}_3$", "Interpreter","latex")
    
    view([30 35])
    
    legend([nanoOrbit, motherOrbit, nanoSat, GMOSat, b_1, negB_1, b_2, b_3], "Nano-satellite Orbit", "Mothercraft Orbit", "Nano-satellite", "Mothercraft", "b_1 axis", "-b_1 axis", "b_2 axis", "b_3 axis", 'Location', 'best')

    % if any(t(k) == tInt)
    %     frames = [frames; getframe(figOrbit)];
    % end
    
    % subplot(1,2,2)
    % % Spacecraft
    % scatter3(out_LMO(k,2), out_LMO(k,3), out_LMO(k,4), 25, 'black', 'filled')
    % scatter3(out_GMO(k,2), out_GMO(k,3), out_GMO(k,4), 25, 'black', 'filled')

    drawnow
    % movieVector = [movieVector; getframe(figOrbit)];
end

% for k = 1:length(frames)
%     switch k
%         case 1
%             imwrite(frames(k).cdata, "pointingAt300s.png")
%         case 2
%             imwrite(frames(k).cdata, "pointingAt2100s.png")
%         case 3
%             imwrite(frames(k).cdata, "pointingAt3400s.png")
%         case 4
%             imwrite(frames(k).cdata, "pointingAt4400s.png")
%         case 5
%             imwrite(frames(k).cdata, "pointingAt5600s.png")
%     end
% end

% movie = VideoWriter('FinalProjectMovie','MPEG-4');
% movie.FrameRate = 30;
% 
% %Open the VideoWriter object, write the movie, and close the file
% open(movie);
% writeVideo(movie, movieVector);
% close(movie);

%% Extract answers to text files
time = t;

t = 300;
sig300 = sig(time == t, :);

t = 2100;
sig2100 = sig(time == t, :);

t = 3400;
sig3400 = sig(time == t, :);

t = 4400;
sig4400 = sig(time == t, :);

t = 5600;
sig5600 = sig(time == t, :);

f1 = fopen("sig300_ans.txt", "w");
ans_sig300 = fprintf(f1, "%.3f %.3f %.3f", sig300(1), sig300(2), sig300(3));
fclose(f1);

f2 = fopen("sig2100_ans.txt", "w");
ans_sig2100 = fprintf(f2, "%.3f %.3f %.3f", sig2100(1), sig2100(2), sig2100(3));
fclose(f2);

f3 = fopen("sig3400_ans.txt", "w");
ans_sig3400 = fprintf(f3, "%.3f %.3f %.3f", sig3400(1), sig3400(2), sig3400(3));
fclose(f3);

f4 = fopen("sig4400_ans.txt", "w");
ans_sig4400 = fprintf(f4, "%.3f %.3f %.3f", sig4400(1), sig4400(2), sig4400(3));
fclose(f4);

f5 = fopen("sig5600_ans.txt", "w");
ans_sig5600 = fprintf(f5, "%.3f %.3f %.3f", sig5600(1), sig5600(2), sig5600(3));
fclose(f5);

