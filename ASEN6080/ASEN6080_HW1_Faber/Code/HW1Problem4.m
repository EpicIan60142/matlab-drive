%% ASEN 6080 HW 1 Problem 4
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
    % Physical parameters
theta0 = deg2rad(122); % Initial spin angle of Earth relative to ECI
elMask = deg2rad(10); % Ground station elevation mask
wEarth = (2*pi)/(24*60*60); % 24 hour Earth spin rate
stations = struct('X0', [], 'Xs', [], 'rho', [], 'rhoDot', [], ...
                  'elAngle', [], 'RU', [], 'fShift', []); % Create struct for storing station info
latitudes = deg2rad([-35.39833; 40.42722; 35.247164]); % [phi_1; phi_2; phi_3]
longitudes = deg2rad([148.981944; 355.749444; 243.205]); % [lambda_1; lambda_2; lambda_3]
colors = ['b', 'r', 'k']; % Station colors for plotting
fT = 8.44e9; % X-band transmission frequency
sigma = 0.5e-6; % range rate noise standard deviation of 0.5 mm/s in km/s

    % Gravitational parameters
mu = 398600; % km^3/s^2
J2 = 1.08264e-3; 
Ri = 6378; % km

    % Orbit parameters (from problem 2)
a = 10000; % km
e = 0.001; 
i = deg2rad(40); % deg -> rad
RAAN = deg2rad(80); % deg -> rad
argPeri = deg2rad(40); % deg -> rad
truAnom = deg2rad(0); % deg -> rad

%     % Orbit parameters (for experimenting, uncomment to apply. 
%     % Note: without changing dt as well, simulations will take drastically longer!)
% a = 42163; 
% e = 1e-5;
% i = rad2deg(1);
% RAAN = deg2rad(70);
% argPeri = deg2rad(0);

    % Z-axis rotation matrix
rotZ = @(theta) [
                    cos(theta), -sin(theta), 0;
                    sin(theta), cos(theta),  0;
                    0,          0,           1
                ];

    % Earth sphere
[earthX_s, earthY_s, earthZ_s] = sphere(20);
earthX_s = Ri*earthX_s;
earthY_s = Ri*earthY_s;
earthX = earthX_s*cos(theta0) - earthY_s*sin(theta0);
earthY = earthX_s*sin(theta0) + earthY_s*cos(theta0);
earthZ = Ri*earthZ_s;
I = imread('2k_earth_daymap.jpg'); % Image from https://www.solarsystemscope.com/textures/

%% Part 4a and 4b
    % Propagate orbit from problem 2

    % Convert orbital elements to cartesian and build initial state
orbital = [mu; a; e; i; RAAN; argPeri; truAnom];
X0_nom = convOrbitalToCartesian(orbital);

    % Part 4e hyperbolic orbit - uncomment to apply
% X0_nom = [-7737.559071593195; -43881.87809094457; 0; 3.347424567061589; 3.828541915617483; 0];

X0_nom = [X0_nom; J2];

    % Propagate for 15 orbits
T = 2*pi*sqrt(a^3/mu);
dt = 10; % sec

tspan = 0:dt:(15*T);
opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

[t_nom, X_nom] = ode45(@(t,X)orbitEOM_MuJ2(t,X,mu,Ri), tspan, X0_nom, opt);

    % Make stations
ax = [];
figure;
hold on; grid on; axis equal;
title("Initial Station States with Spacecraft Orbit")
for k = 1:length(latitudes)
    spherical = [Ri; latitudes(k); longitudes(k)];
    body = toBodyFromSpherical(spherical);
    stations(k).X0 = rotZ(theta0)*body;
    ax(k) = plot3(stations(k).X0(1), stations(k).X0(2), stations(k).X0(3), 'Color', colors(k), 'Marker', '.', 'MarkerSize', 15);
end
earth = surf(earthX, earthY, earthZ, 'FaceAlpha', 0.5);
set(earth,'FaceColor','texturemap','cdata',I,'edgecolor','none');
view([30 35])
ax(k+1) = plot3(X_nom(:,1), X_nom(:,2), X_nom(:,3), 'm-');
legend(ax, ["Station 1", "Station 2", "Station 3", "Spacecraft Orbit"])

    % Propagate station states
for k = 1:length(t_nom)
    dTheta = wEarth*t_nom(k);
    for kk = 1:length(stations)
        r = rotZ(dTheta)*stations(kk).X0;
        v = cross([0;0;wEarth], rotZ(dTheta)*stations(kk).X0);

        stations(kk).Xs = [stations(kk).Xs; [r', v']];
        
        y = generateRngRngRate(X_nom(k,:), stations(kk).Xs(k,:), elMask);
        
        stations(kk).rho = [stations(kk).rho; y(1)];
        stations(kk).rhoDot = [stations(kk).rhoDot; y(2)];
        stations(kk).elAngle = [stations(kk).elAngle; y(3)];
    end
end

    % Create measurement plots
figure; tl = tiledlayout(3,1); ax = [];
title(tl, "Ground Station Measurements vs. Time - \rho and \rho-Rate")
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        ax(k) = plot(t_nom, stations(k).rho, 'Color', colors(k), 'Marker', 'x', 'MarkerSize', 2);
    end
    xlabel("Time [sec]"); ylabel("\rho [km]")
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        plot(t_nom, stations(k).rhoDot, 'Color', colors(k), 'Marker', 'o', 'MarkerSize', 2);
    end
    xlabel("Time [sec]"); ylabel("\rhoDot [km/s]")
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        plot(t_nom, rad2deg(stations(k).elAngle), 'Color', colors(k), 'Marker', 'square', 'MarkerSize', 2);
    end
    yline(rad2deg(elMask),'k--')
    xlabel("Time [sec]"); ylabel("Elevation Angle [deg]")
legend(ax, ["Station 1", "Station 2", "Station 3"], 'location', 'bestoutside')

%% Part 4c
    % Retake measurements with varying fT
for k = 1:length(t_nom)
    dTheta = wEarth*t_nom(k);
    for kk = 1:length(stations)
        
        y = generateRUDoppler(X_nom(k,:), stations(kk).Xs(k,:), elMask, fT*kk);
        
        stations(kk).RU = [stations(kk).RU; y(1)];
        stations(kk).fShift = [stations(kk).fShift; y(2)];
    end
end

    % Plot new measurements
figure; tl = tiledlayout(3,1); ax = [];
title(tl, "Ground Station Measurements vs. Time - RU and f_{shift}")
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        ax(k) = plot(t_nom, stations(k).RU, 'Color', colors(k), 'Marker', 'x', 'MarkerSize', 2);
    end
    xlabel("Time [sec]"); ylabel("RU")
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        plot(t_nom, stations(k).fShift, 'Color', colors(k), 'Marker', 'o', 'MarkerSize', 2);
    end
    xlabel("Time [sec]"); ylabel("f_{shift} [Hz]")
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        plot(t_nom, rad2deg(stations(k).elAngle), 'Color', colors(k), 'Marker', 'square', 'MarkerSize', 2);
    end
    yline(rad2deg(elMask),'k--')
    xlabel("Time [sec]"); ylabel("Elevation Angle [deg]")
legend(ax, ["Station 1 (8.44 GHz)", "Station 2 (16.88 GHz)", "Station 3 (25.32 GHz)"], 'location', 'bestoutside')

%% Part 4d

    % Add noise to range-rate measurements
for kk = 1:length(stations)
    noise = mvnrnd(0, sigma, length(t_nom)); % Generate noise
    stations(kk).rhoDotNoisy = stations(kk).rhoDot + noise;
end

    % Create measurement plots
figure; tl = tiledlayout(3,1); ax = [];
title(tl, "Ground Station Measurements vs. Time - Noisy \rho-Rate")
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        ax(k) = plot(t_nom, stations(k).rho, 'Color', colors(k), 'Marker', 'x', 'MarkerSize', 2);
    end
    xlabel("Time [sec]"); ylabel("\rho [km]")
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        plot(t_nom, stations(k).rhoDotNoisy, 'Color', colors(k), 'Marker', 'o', 'MarkerSize', 2);
    end
    xlabel("Time [sec]"); ylabel("\rhoDot [km/s]")
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        plot(t_nom, rad2deg(stations(k).elAngle), 'Color', colors(k), 'Marker', 'square', 'MarkerSize', 2);
    end
    yline(rad2deg(elMask),'k--')
    xlabel("Time [sec]"); ylabel("Elevation Angle [deg]")
legend(ax, ["Station 1", "Station 2", "Station 3"], 'location', 'bestoutside')

    % Plot difference in original and noisy data
figure; tl = tiledlayout(3,1); ax = [];
title(tl, "Noisy minus Original Ground Station Measurements vs. Time")
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        ax(k) = plot(t_nom, stations(k).rho - stations(k).rho, 'Color', colors(k), 'Marker', 'x', 'MarkerSize', 2);
    end
    xlabel("Time [sec]"); ylabel("\Delta\rho [km]")
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        plot(t_nom, stations(k).rhoDotNoisy - stations(k).rhoDot, 'Color', colors(k), 'Marker', 'o', 'MarkerSize', 2, 'LineStyle', 'none');
    end
    xlabel("Time [sec]"); ylabel("\Delta\rhoDot [km/s]")
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        plot(t_nom, rad2deg(stations(k).elAngle - stations(k).elAngle), 'Color', colors(k), 'Marker', 'square', 'MarkerSize', 2);
    end
    xlabel("Time [sec]"); ylabel("\Delta Elevation Angle [deg]")
legend(ax, ["Station 1", "Station 2", "Station 3"], 'location', 'bestoutside')

%% Fun!
return; % Comment out to animate!

    % Animate 4a-b!
movieVector = [];
dTime = 10; % How many frames to skip to speed up simulation
fig = figure('Position',[0 0 1920 1080]); 
for k = 1:dTime:length(t_nom)
    clf; 

    tl = tiledlayout(3,2);
    title(tl, "ASEN 6080 HW 1 Animation")

        % 3D view
    nexttile([3 1]) % Make 3d plot span 3 rows and 1 column
        hold on; grid on; axis equal;
    
        titleText = sprintf("Scenario at t = %.0f sec", t_nom(k));
        title(titleText);
    
            % Satellite orbit
        plot3(X_nom(:,1), X_nom(:,2), X_nom(:,3), 'c--');
    
        visible = false;
        for kk = 1:length(stations)
            if ~isnan(stations(kk).elAngle(k))
                visible = true;
            end
        end
        
        if visible
            plot3(X_nom(k,1), X_nom(k,2), X_nom(k,3), 'm.', 'MarkerSize', 25);
        else
            plot3(X_nom(k,1), X_nom(k,2), X_nom(k,3), 'k.', 'MarkerSize', 15);
        end
    
            % Ground stations
        for kk = 1:length(stations)
            plot3(stations(kk).Xs(k,1), stations(kk).Xs(k,2), stations(kk).Xs(k,3), 'Color', colors(kk), 'Marker', '.', 'MarkerSize', 15);
        end
    
            % Earth
        earthX = earthX_s*cos(theta0 + wEarth*t_nom(k)) - earthY_s*sin(theta0 + wEarth*t_nom(k));
        earthY = earthX_s*sin(theta0 + wEarth*t_nom(k)) + earthY_s*cos(theta0 + wEarth*t_nom(k));
        earth = surf(earthX, earthY, earthZ, 'FaceAlpha', 0.5);
        set(earth,'FaceColor','texturemap','cdata',I,'edgecolor','none');
    
            % Format figure
        view([30 35])
        xlim([-1.1*a, 1.1*a]);
        ylim([-1.1*a, 1.1*a]);
        zlim([-1.1*a, 1.1*a]);

        % Measurements
    nexttile
        hold on; grid on;
        title("Ground Station measurements")
        for kk = 1:length(stations)
            ax(kk) = plot(t_nom, stations(kk).rho, 'Color', colors(kk), 'Marker', 'x', 'MarkerSize', 2);
            plot(t_nom(k), stations(kk).rho(k), 'Color', colors(kk), 'Marker', '.', 'MarkerSize', 15);
        end
        xlabel("Time [sec]"); ylabel("\rho [km]")
    nexttile
        hold on; grid on;
        for kk = 1:length(stations)
            plot(t_nom, stations(kk).rhoDot, 'Color', colors(kk), 'Marker', 'o', 'MarkerSize', 2);
            plot(t_nom(k), stations(kk).rhoDot(k), 'Color', colors(kk), 'Marker', '.', 'MarkerSize', 15);
        end
        xlabel("Time [sec]"); ylabel("\rhoDot [km/s]")
    nexttile
        hold on; grid on;
        for kk = 1:length(stations)
            plot(t_nom, rad2deg(stations(kk).elAngle), 'Color', colors(kk), 'Marker', 'square', 'MarkerSize', 2);
            plot(t_nom(k), rad2deg(stations(kk).elAngle(k)), 'Color', colors(kk), 'Marker', '.', 'MarkerSize', 15);
        end
        yline(rad2deg(elMask),'k--')
        xlabel("Time [sec]"); ylabel("Elevation Angle [deg]")
    legend(ax, ["Station 1", "Station 2", "Station 3"], 'location', 'bestoutside')

        % Update figure and make movie vector
    drawnow
    movieVector = [movieVector; getframe(fig)];

end

%     % Make movie
% movie = VideoWriter('ASEN6080_HW1_New','MPEG-4');
% movie.FrameRate = 30;
% 
%     % Open the VideoWriter object, write the movie, and close the file
% open(movie);
% writeVideo(movie, movieVector);
% close(movie);
