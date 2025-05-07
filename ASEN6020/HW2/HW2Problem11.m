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
dTheta = linspace(0.001,2*pi-0.001,10); % rad
TOF = linspace(0.1,4*24*3600,200); % sec, 0 - 4 day transfer time

%% Solve Lambert's Problem for every r, dTheta combo

fprintf("\nFinding transfers for solution space\n")

if false
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
    save("HW2Problem11Data-new.mat",'transfers','-mat');
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
NO = @(theta, inc, RAAN)EA2DCM([-theta, -inc, -RAAN], [3, 1, 3]); % orbital -> inertial

    % Plot all transfers
for ratio = r
    figure; hold on; grid on; axis equal; colors = 'jet';
    titleText = sprintf("Elliptical Transfers from r_1 = %.3e km to r_2 = %.3e km \n given \\Delta\\theta* \\epsilon [%.3f, %.3f] and r \\epsilon [%.3f, %.3f]", r1, ratio*r1, dTheta(1), dTheta(end), r(1), r(end));
    title(titleText)
    for k = 1:length(transfers)
        if abs(transfers(k).diagnostic.r2 - ratio*r1) < 1e-10 % Only plot transfers that don't hit Earth
            a = transfers(k).a;
            e = transfers(k).e;
            inc = transfers(k).orbElems.inc;
            RAAN = transfers(k).orbElems.RAAN;
            argPeri = transfers(k).orbElems.argPeri;
            if transfers(k).type == "Elliptical"
                TA = deg2rad(real(transfers(k).TA1_deg):0.5:real(transfers(k).TA2_deg));
                theta = TA + argPeri;
            else
                continue
                % a = transfers(534).a;
                % e = transfers(534).e;
                % inc = transfers(534).orbElems.inc;
                % RAAN = transfers(534).orbElems.RAAN;
                % argPeri = transfers(534).orbElems.argPeri;
                % TA1 = deg2rad(transfers(534).TA1_deg); %wrapTo2Pi(deg2rad(transfers(534).TA1_deg));
                % TA2 = deg2rad(transfers(534).TA2_deg); %wrapTo2Pi(deg2rad(transfers(534).TA2_deg));
                % TA = TA1:deg2rad(0.5):TA2;
                % theta = wrapToPi(TA + argPeri);
            end
            
            radius = (a*(1-e^2))./(1+e*cos(TA));
            if any(radius < 6378) % Don't plot if the transfer intersects Earth
                continue
            end
                % rhat-thetahat-hhat frame
            R = [radius;zeros(size(radius));zeros(size(radius))];
                % rhat-thetahat-hhat -> Cartesian frame
            for kk = 1:length(theta)
                R(:,kk) = real(NO(theta(kk),inc,RAAN)*R(:,kk));
            end
            % plot3(R(1,:), R(2,:), R(3,:))
            % if transfers(k).dV_mag_total% < 80 % Limit plotting based on delta V's
                color_line3d(real(transfers(k).TOF_calc)*ones(size(R(1,:))), R(1,:), R(2,:), R(3,:)); % Create colored line based on total delta V required
                startTrans = plot3(R(1,1), R(2,1), R(3,1), 'g.', 'markerSize', 15); % Start of transfer
                stopTrans = plot3(R(1,end), R(2,end), R(3,end), 'r.', 'markerSize', 15); % End of transfer
            % end
        end
    end

        % Plot original orbits
    for k = 1:2
        TA = 0:0.01:2*pi;
        switch k
            case 1
                e = 0;
                a = r1;
                radius = (a*(1-e^2))./(1+e*cos(TA));
                theta = TA;
                    % rhat-thetahat-hhat frame
                R = [radius;zeros(size(radius));zeros(size(radius))];
                    % rhat-thetahat-hhat -> Cartesian frame
                for kk = 1:length(theta)
                    R(:,kk) = real(NO(theta(kk),inc,RAAN)*R(:,kk));
                end
            case 2
                e = 0;
                a = ratio*r1;
                radius = (a*(1-e^2))./(1+e*cos(TA));
                theta = TA;
                    % rhat-thetahat-hhat frame
                R = [radius;zeros(size(radius));zeros(size(radius))];
                    % rhat-thetahat-hhat -> Cartesian frame
                for kk = 1:length(theta)
                    R(:,kk) = real(NO(theta(kk),inc,RAAN)*R(:,kk));
                end
        end
        plot3(R(1,:), R(2,:), R(3,:), 'k--')
    end
    
        % Earth
    e = 0;
    a = 6378;
    radius = (a*(1-e^2))./(1+e*cos(TA));
    theta = TA;
        % rhat-thetahat-hhat frame
    R = [radius;zeros(size(radius));zeros(size(radius))];
        % rhat-thetahat-hhat -> Cartesian frame
    for kk = 1:length(theta)
        R(:,kk) = real(NO(theta(kk),inc,RAAN)*R(:,kk));
    end
    plot3(R(1,:), R(2,:), R(3,:), 'k-')

    try
        c = colorbar; c.Label.String = "Time of Flight [sec]"; colormap(colors)
        xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
        xlim([-3e5, 3e5]), ylim([-3e5, 3e5]); %zlim([-5e6, 5e6]); %view([30 35])
        legend([startTrans, stopTrans], ["Start of transfer", "End of transfer"], "location", 'best')
        drawnow
    catch
        fprintf("\nTransfers couldn't match TOF to machine precision\n")
    end
end





