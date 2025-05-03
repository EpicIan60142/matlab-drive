%% ASEN 6020 Final Project Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Logistics
addpath(genpath(".\CubeSats"));
addpath(genpath(".\RaceCourse"));
addpath(genpath(".\Plotting"));

% %% Enable parallel computing
% parpool;

%% Ask user for input
menuMsg = sprintf("\n\t\tWelcome to the CubeSat racer! What would you like to do?\n\n" + ...
                  "\t1. Generate and run a new race course (time intensive!)\n" + ...
                  "\t2. Analyze an old run\n" + ...
                  "\t3. Exit\n\n" + ...
                  "\tChoice (enter number of option to do): ");
choice = input(menuMsg);

%% Decide what to do
switch choice
    case 1
        %% Setup Problem Pameters
        fprintf("Generating and solving new race course, this will take some time!\n")

            % Define ring semi-major and minor axis ranges
        semiMaj = [1, 5]; % m
        semiMin = [1, 5]; % m
        
            % Define inter-ring parameter ranges
        interRingDist = [300, 350]; % m
        azimuthAngle = deg2rad([-90, 90]); % deg -> rad
        elevationAngle = deg2rad([-90, 90]); % deg -> rad
        
            % Define number of rings range
        numRings = [15, 25];
        
            % Define course origin orbit
        T = 180*60; % 180 minute circular orbit
        n = 2*pi/T;
        
            % Define course parameters structure
        courseParams = struct("semiMaj", semiMaj, "semiMin", semiMin, "dist", interRingDist, ...
                              "azAng", azimuthAngle, "elAng", elevationAngle, "numRings", numRings, ...
                              "n", n);
        
            % ODE45 options
        opt = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
        
        %% Generate race course
            % Set rng seed for testing
        if true
            seed = 2*69420; % cool options: 0, 2, 3, 4
            rng(seed);
        end
        
            % Make course
        rings = generateRaceCourse(courseParams);
        
            % Define CubeSat starting ring as the maximum distance behind the start
            % ring at a 3 sigma covariance of 5x the largest semi-major/minor axis
        startRing = generateRing(50*max(semiMaj), 50*max(semiMin), 0, 0, -max(interRingDist), rings(1));
        rings(1).params.lastRing = startRing;

            % Define CubeSat ending ring as the minimum distance after the last 
            % ring at a covariance corresponding to the smallest possible ring
        endRing = generateRing(min(semiMaj),min(semiMin),0,0,min(interRingDist),rings(end));
        rings = [rings; endRing];
        
            % Move race course so origin is the mean of all ring centers
        centers = [];
        for k = 1:length(rings)
            centers = [centers, rings.center];
        end
        
        newOrigin = mean(centers, 2);
        
        startRing.center = startRing.center - newOrigin;
        startRing.params.lastRing.center = startRing.params.lastRing.center - newOrigin;
        for k = 1:length(rings)
            rings(k).center = rings(k).center - newOrigin;
            rings(k).params.lastRing.center = rings(k).params.lastRing.center - newOrigin;
        end
        endRing.center = endRing.center - newOrigin;
        endRing.params.lastRing.center = endRing.params.lastRing.center - newOrigin;
        
            % Plot race course
        figure(420)
        hold on; grid on; axis equal
        title("Generated Race Course")
        for k = 1:length(rings)-1
            scatter3(rings(k).center(1), rings(k).center(2), rings(k).center(3), 10, k, 'filled')
            quiver3(rings(k).center(1), rings(k).center(2), rings(k).center(3), rings(k).normal(1), rings(k).normal(2), rings(k).normal(3), 10, 'filled', 'k-')
            plotRing(rings(k), 'k-');
        end
        cubeStart = plotRing(startRing, 'g-');
        quiver3(startRing.center(1), startRing.center(2), startRing.center(3), startRing.normal(1), startRing.normal(2), startRing.normal(3), 10, 'filled', 'k-')
        cubeEnd = plotRing(endRing, 'r-');
        
        courseCenter = scatter3(0, 0, 0, 20, 'k', 'filled', 'h');
        
        xlabel("X [m]"); ylabel("Y [m]"); zlabel("Z [m]"); cBar = colorbar;
        cBar.Label.String = "Ring Number"; colormap("cool");
        
        %% Generate CubeSats
            % Max thrust based on https://www.cubesatshop.com/product/cluster-ifm-nano-thruster-smallsats/ - 10-500 uN
        numSats = 4;
        names = ["Eeny", "Meeny", "Miny", "Mo"]; % Names :P Mo has attitude control system problems, and can only use x-y-z thrusters!
        markers = ["o", "^", "square", "diamond"]; % Markers for plotting
        colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"]; % Colors for plotting
        thrusts = [100, 250, 500, 500]*(10^-3); % Maximum thrust values
        thrustConfigs = [false, false, false, true]; % Whether thrust is split amongst x-y-z thrusters or a general direction
        
            % Make CubeSats

        cubesats = [];
        for k = 1:numSats
            cubesats = [cubesats; generateCubesat(thrusts(k), thrustConfigs(k), startRing, names(k), markers(k), colors(k))];
        end
        
            % Plot starting positions
        cubeAx = []; cubeLabels = [];
        for k = 1:length(cubesats)
            cubeAx = [cubeAx, plot3(cubesats(k).X0(1), cubesats(k).X0(2), cubesats(k).X0(3), 'Color', cubesats(k).color, 'Marker', cubesats(k).marker, 'MarkerFaceColor', 'k')];
            cubeLabels = [cubeLabels, sprintf("CubeSat %s", cubesats(k).name)];
        end
        
        legend([cubeStart, cubeEnd, courseCenter, cubeAx], ["CubeSat 3\sigma Starting Sample Space", "CubeSat End Target Space", "Race Course Origin", cubeLabels], 'location', 'best');
        
            % Show off course
        for k = 0:360
            view(-30 + k, 35);
            drawnow
        end
        
        % Test dynamics
        % X0 = [cubesats(1).X0; 1e-3*ones(6,1)];
        % [t,XTest] = ode45(@(t,X)CHWEOM(t,X,cubesats(1),courseParams), [0 T], X0, opt);
        
        % figure
        % hold on; grid on;
        % plot3(XTest(:,1), XTest(:,2), XTest(:,3));
        % view([30 35])
        
        % return;
        
        %% Run the race course!
        debug = [false; false; false]; % iteration plotting+display; trajectory plotting; disable progress sequence
        for k = 1:length(cubesats)
            fprintf("\nCubesat %s is starting the race course!\n", cubesats(k).name)
            for kk = 1:length(rings)-1 % Intermediate ring problem
                    % Solve trajectory
                [t,X,optParams,initGuess,cost] = solveTrajectory_int(cubesats(k), rings(kk), courseParams, opt, debug);
                cubesats(k).X = [cubesats(k).X; X];
                cubesats(k).t = [cubesats(k).t; t];
                cubesats(k).optParams = [cubesats(k).optParams, optParams]; % Parameters solved by fmincon
                cubesats(k).initGuess = [cubesats(k).initGuess, initGuess]; % Initial guess for parameters
                cubesats(k).ringSeg = [cubesats(k).ringSeg, [kk-1; kk]]; % Ring segment of this trajectory (from kk-1 to kk)
                cubesats(k).tSeg = [cubesats(k).tSeg, [t(1); t(end)]]; % Time this ring segment started and ended
                cubesats(k).cost = [cubesats(k).cost, cost]; % Cost for this ring segment

                    % Update t0 and X0
                cubesats(k).X0 = X(end,1:6)';
                cubesats(k).t0 = t(end);
        
                    % Plot segment
                titleText = sprintf("Cubesat trajectory segment: Ring %.0f to %.0f", kk-1, kk);
                xLabel = "Radial [km]"; yLabel = "Along Track [km]"; zLabel = "Cross Track [km]";
                plotSegment(cubesats(k), rings(kk), t, X, kk + 1, titleText, xLabel, yLabel, zLabel);
                
                    % Report progress
                fprintf("\tCubesat %s passed through ring %.0f in %.3f sec!\n", cubesats(k).name, kk, t(end) - t(1));
            end
            fprintf("\n\tCubesat %s finished the course in %.3f sec!\n", cubesats(k).name, cubesats(k).t(end)-cubesats(k).t(1));
                % Add trajectory to race plot
            figure(420)
            % title(sprintf("Cubesat %s Race Course Trajectory", cubesats(k).name));
            plot3(cubesats(k).X(:,1), cubesats(k).X(:,2), cubesats(k).X(:,3), '-', 'Color', cubesats(k).color, 'DisplayName', sprintf("Cubesat %s trajectory", cubesats(k).name));
        
            figure; tl = tiledlayout(1,2);
            title(tl, sprintf("Cubesat %s Adjoints", cubesats(k).name))
            nexttile
                title(sprintf("Cubesat %s Position Adjoints", cubesats(k).name))
                hold on; grid on; axis equal;
                plot3(cubesats(k).X(:,7), cubesats(k).X(:,8), cubesats(k).X(:,9), 'k.');
                view([30 35])
            nexttile
                title(sprintf("Cubesat %s Velocity Adjoints", cubesats(k).name))
                hold on; grid on; axis equal;
                plot3(cubesats(k).X(:,10), cubesats(k).X(:,11), cubesats(k).X(:,12), 'k.');
                view([30 35])
        end
        
            % Save run!
        time = datetime('now');
        runName = sprintf("Run_%s", time);
        runName = replace(runName, " ", "_");
        runName = replace(runName, ":", "");
        runName = sprintf(".\\Runs\\%s.mat", runName);
        save(runName, "cubesats", "rings", "startRing", "endRing", "courseParams", '-mat');

    case 2
        %% Choose run to analyze
        fprintf("\n\tChoose a run to analyze!\n");
        cd .\Runs\
        [file, path, ~] = uigetfile('.mat', "Pick a MAT file to analyze");
        cd ..\

        load([path,file])

            % Plot race course
        figure(420)
        hold on; grid on; axis equal
        title("Generated Race Course")
        for k = 1:length(rings)-1
            scatter3(rings(k).center(1), rings(k).center(2), rings(k).center(3), 10, k, 'filled')
            quiver3(rings(k).center(1), rings(k).center(2), rings(k).center(3), rings(k).normal(1), rings(k).normal(2), rings(k).normal(3), 10, 'filled', 'k-')
            plotRing(rings(k), 'k-');
        end
        cubeStart = plotRing(startRing, 'g-');
        quiver3(startRing.center(1), startRing.center(2), startRing.center(3), startRing.normal(1), startRing.normal(2), startRing.normal(3), 10, 'filled', 'k-')
        cubeEnd = plotRing(endRing, 'r-');
        
        courseCenter = scatter3(0, 0, 0, 20, 'k', 'filled', 'h');
        
        xlabel("X [m]"); ylabel("Y [m]"); zlabel("Z [m]"); cBar = colorbar;
        cBar.Label.String = "Ring Number"; colormap("cool"); view([30 35]);

            % Plot starting positions
        cubeAx = []; cubeLabels = [];
        for k = 1:length(cubesats)
            cubeAx = [cubeAx, plot3(cubesats(k).X(1,1), cubesats(k).X(1,2), cubesats(k).X(1,3), 'Color', cubesats(k).color, 'Marker', cubesats(k).marker, 'MarkerFaceColor', 'k')];
            cubeLabels = [cubeLabels, sprintf("CubeSat %s", cubesats(k).name)];
        end
        
        legend([cubeStart, cubeEnd, courseCenter, cubeAx], ["CubeSat 3\sigma Starting Sample Space", "CubeSat End Target Space", "Race Course Origin", cubeLabels], 'location', 'best');

            % Plot trajectory segments
        for k = 1:length(cubesats)
            fprintf("\nPlotting Cubesat %s's Trajectory!\n", cubesats(k).name)
            for kk = 1:length(rings)-1 % Intermediate ring problem
                    % Find indices
                kStart = find(cubesats(k).t == cubesats(k).tSeg(1,kk), 1, 'first');
                kEnd = find(cubesats(k).t == cubesats(k).tSeg(2,kk), 1, 'last');

                    % Plot segment
                figure(kk)
                title(sprintf("Cubesat trajectory segment: Ring %.0f to %.0f", kk-1, kk));
                hold on; grid on; axis equal
                if kk == 1
                    plotRing(startRing, 'g-');
                else
                    plotRing(rings(kk).params.lastRing, 'g-');
                end
                plot3(cubesats(k).X(kStart,1), cubesats(k).X(kStart,2), cubesats(k).X(kStart,3), 'k', 'Marker', cubesats(k).marker);
                plotRing(rings(kk), 'r-');
                plot3(cubesats(k).X(kStart:kEnd,1), cubesats(k).X(kStart:kEnd,2), cubesats(k).X(kStart:kEnd,3), '--', 'Color', cubesats(k).color);
                view([30 35]);
        
                    % Report progress
                fprintf("\tCubesat %s passed through ring %.0f in %.3f sec!\n", cubesats(k).name, kk, cubesats(k).t(kEnd) - cubesats(k).t(kStart));
            end
            fprintf("\n\tCubesat %s finished the course in %.3f sec!\n", cubesats(k).name, cubesats(k).t(end)-cubesats(k).t(1));
            
                % Add trajectory to race course plot
            figure(420)
            plot3(cubesats(k).X(:,1), cubesats(k).X(:,2), cubesats(k).X(:,3), '-', 'Color', cubesats(k).color, 'DisplayName', sprintf("Cubesat %s trajectory", cubesats(k).name));
        
            figure; tl = tiledlayout(1,2);
            title(tl, sprintf("Cubesat %s Adjoints", cubesats(k).name))
            nexttile
                title(sprintf("Cubesat %s Position Adjoints", cubesats(k).name))
                hold on; grid on; axis equal;
                plot3(cubesats(k).X(:,7), cubesats(k).X(:,8), cubesats(k).X(:,9), 'k.');
                view([30 35])
            nexttile
                title(sprintf("Cubesat %s Velocity Adjoints", cubesats(k).name))
                hold on; grid on; axis equal;
                plot3(cubesats(k).X(:,10), cubesats(k).X(:,11), cubesats(k).X(:,12), 'k.');
                view([30 35])
        end

    case 3
        %% Say goodbye
        fprintf("\n\tHave a great day!\n")
        return;
end

