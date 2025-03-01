%% ASEN 6020 HW 2 Problem 11 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Define constants and solution space
    % Constants
r1 = 10000; % km
mu = 398600.4415; % km^3/s^2
v_lc1 = sqrt(mu/r1);
X_1 = [r1; 0; 0; 0; v_lc1; 0]; % Initial state is a constant with theta_1 = 0!

    % Solution space
r = linspace(11.94,25,10);
dTheta = linspace(0.001,2*pi,10); % rad
TOF = linspace(0.1,24*3600,100); % sec, 0 - 1 day transfer time

%% Solve Lambert's Problem for every r, dTheta combo

fprintf("\nFinding transfers for solution space\n")

if true
    transfers = [];
    for t = TOF
        for k = 1:length(r)
            for kk = 1:length(dTheta)
               fprintf("\nTOF = %.3e sec, r = %.3f, dTheta = %.3f rad, k = %.0f, kk = %.0f\n", t, r(k), dTheta(kk), k, kk);
    
                    % Pull out parameters for this run
                r_kk = r(k);
                dTheta_kk = dTheta(kk);
                    
                    % Determine final parameters
                r2 = r_kk*r1;
                v_lc2 = sqrt(mu/(r2));
        
                    % Rotate from r,theta,h frame to x,y,z frame
                X_2 = blkdiag(rotZ(dTheta_kk), rotZ(dTheta_kk))*[r2; 0; 0; 0; v_lc2; 0];
        
                    % Determine if transfer arc is long or short
                if dTheta_kk > pi
                    lt180 = 0;
                else
                    lt180 = 1;
                end
    
                    % Find optimal transfer that connects X_1 and X_2
                transfers = [transfers; solveLambertsProblem(X_1, X_2, t, lt180, mu)];
        
            end
        end
    end
    save("HW2Problem11Data.mat",'transfers','-mat');
else
    load("HW2Problem11Data.mat")
end

%% Plot Pareto Front
fprintf("\nPlotting Pareto Front\n")

figure;
semilogy(transfers(1).TOF_calc, transfers(1).dV_mag_total)
hold on; grid on;
title("Pareto Front between TOF and \DeltaV")
for k = 2:length(transfers)
    semilogy(real(transfers(k).TOF_calc), transfers(k).dV_mag_total, 'b.')
end
xlabel("Time of Flight [sec]"); ylabel("\DeltaV [km/s]")
drawnow;

% return

%% Plot transfers
fprintf("\nPlotting Transfers\n")
    % DCM function handle
% NO = @(theta, inc, RAAN)EA2DCM([-theta, -inc, -RAAN], [3, 1, 3]); % orbital -> inertial

    % Plot all transfers
for T = TOF
    figure; hold on; grid on; colors = 'jet';
    titleText = sprintf("All Possible Transfers from r_1 = %.3e km to varying r_2, \n given TOF = %.3e sec", r1, rad2deg(T));
    title(titleText)
    for k = 1:length(transfers)
        if abs(transfers(k).TOF_calc - T) < 1e-10
            a = transfers(k).a;
            e = transfers(k).e;
            % inc = transfers(k).orbElems.inc;
            % RAAN = transfers(k).orbElems.RAAN;
            % argPeri = transfers(k).orbElems.argPeri;
            TA = deg2rad(real(transfers(k).TA1_deg):0.5:real(transfers(k).TA2_deg));
            theta = TA;% + argPeri;
            r = (a*(1-e^2))./(1+e*cos(TA));
            % rhat-thetahat-hhat frame
            R = [r;zeros(size(r));zeros(size(r))];
            % rhat-thetahat-hhat -> Cartesian frame
            for kk = 1:length(theta)
                R(:,kk) = real(rotZ(theta(kk))*R(:,kk));
            end
            % plot3(R(1,:), R(2,:), R(3,:))
            if transfers(k).dV_mag_total% < 80 % Limit plotting based on delta V's
                color_line3d(transfers(k).dV_mag_total*ones(size(R(1,:))), R(1,:), R(2,:), R(3,:)); % Create colored line based on total delta V required
                startTrans = plot3(R(1,1), R(2,1), R(3,1), 'g.', 'markerSize', 15); % Start of transfer
                stopTrans = plot3(R(1,end), R(2,end), R(3,end), 'r.', 'markerSize', 15); % End of transfer
            end
        end
    end
    c = colorbar; c.Label.String = "Total Delta V [km/s]"; colormap(colors)
    xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
    xlim([-5e6, 5e6]), ylim([-5e6, 5e6]); zlim([-5e6, 5e6]); %view([30 35])
    legend([startTrans, stopTrans], ["Start of transfer", "End of transfer"], "location", 'best')
    drawnow
end





