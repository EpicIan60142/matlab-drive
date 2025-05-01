function animateProblem(X, stations, pConst, titleTextTotal, movieTitle, saveMovie)
% Function that animates the results of an OD problem.
%   Inputs:
%       - X: Full state history of the problem, containing the spacecraft 
%            state at each instant in time. Organized as
%            [X_1, X_2, ..., X_tMeas]
%       - stations: Stations struct as defined by makeStations.m, must be
%                   populated with measurements and measurement times
%       - titleTextTotal: Overall title for the animation figure
%       - movieTitle: Name of the video to save the animation to (.mp4)
%       - saveMovie: Boolean indicating whether or not to save the movie to
%                    a file
%   Outputs:
%       - None
%
%   By: Ian Faber, 02/15/2025
%

%% Set up scenario
    % Get unified measurement time vector
[tMeas, ~, ~, Xs, ~] = processStations(stations);
tMeas = tMeas(2:end);
Xs = Xs(2:end);

    % Extract states to plot
X_sc = X;
X_stat = cellfun(@transpose,Xs,'UniformOutput',false);

    % Animation setup
movieVector = [];
dTime = 1; % How many frames to skip to speed up simulation
theta0 = deg2rad(122);

    % Earth sphere
Ri = pConst.Ri;
wEarth = pConst.wPlanet;
[earthX_s, earthY_s, earthZ_s] = sphere(20);
earthX_s = Ri*earthX_s;
earthY_s = Ri*earthY_s;
earthX = earthX_s*cos(theta0) - earthY_s*sin(theta0);
earthY = earthX_s*sin(theta0) + earthY_s*cos(theta0);
earthZ = Ri*earthZ_s;
I = imread('2k_earth_daymap.jpg'); % Image from https://www.solarsystemscope.com/textures/

%% Animate!
fig = figure('Position',[0 0 1920 1080]); 
for k = 1:dTime:length(tMeas)
    clf; 

    tl = tiledlayout(2,2);
    title(tl, titleTextTotal)

        % 3D view
    nexttile([2 1]) % Make 3d plot span 2 rows and 1 column
        hold on; grid on; axis equal;
    
        titleText = sprintf("Scenario at t = %.0f sec", tMeas(k));
        title(titleText);
    
            % Satellite orbit
        plot3(X_sc(1,:), X_sc(2,:), X_sc(3,:), 'c.');
        
            % Satellite at time of measurement
        plot3(X_sc(1,k), X_sc(2,k), X_sc(3,k), 'm.', 'MarkerSize', 25);
    
            % Ground stations
        for kk = 1:1%length(stations)
            stat = X_stat{k};
            plot3(stat(3*kk-2), stat(3*kk-1), stat(3*kk), 'Color', stations(kk).color, 'Marker', '.', 'MarkerSize', 15);
        end
    
            % Earth
        earthX = earthX_s*cos(theta0 + wEarth*tMeas(k)) - earthY_s*sin(theta0 + wEarth*tMeas(k));
        earthY = earthX_s*sin(theta0 + wEarth*tMeas(k)) + earthY_s*cos(theta0 + wEarth*tMeas(k));
        earth = surf(earthX, earthY, earthZ, 'FaceAlpha', 0.5);
        set(earth,'FaceColor','texturemap','cdata',I,'edgecolor','none');
    
            % Format figure
        % view([30+k/3, 35])
        view([30, 35])
        xlabel("X [m]"); ylabel("Y [m]"); zlabel("Z [m]");
        % xlim([-1.1*a, 1.1*a]);
        % ylim([-1.1*a, 1.1*a]);
        % zlim([-1.1*a, 1.1*a]);

        % Measurements
    nexttile
        hold on; grid on;
        title("Ground Station measurements")
        labels = [];
        for kk = 1:length(stations)
            ax(kk) = plot(stations(kk).tMeas, stations(kk).rho, 'Color', stations(kk).color, 'Marker', 'x', 'MarkerSize', 2, 'LineStyle', 'none');
            mask = stations(kk).tMeas == tMeas(k);
            plot(stations(kk).tMeas(mask), stations(kk).rho(mask), 'Color', stations(kk).color, 'Marker', '.', 'MarkerSize', 15);
            labels = [labels; sprintf("Station %.0f", stations(kk).id)];
        end
        xlabel("Time [sec]"); ylabel("\rho [m]")
    nexttile
        hold on; grid on;
        for kk = 1:length(stations)
            plot(stations(kk).tMeas, stations(kk).rhoDot, 'Color', stations(kk).color, 'Marker', 'o', 'MarkerSize', 2, 'LineStyle', 'none');
            mask = stations(kk).tMeas == tMeas(k);
            plot(stations(kk).tMeas(mask), stations(kk).rhoDot(mask), 'Color', stations(kk).color, 'Marker', '.', 'MarkerSize', 15);
        end
        xlabel("Time [sec]"); ylabel("\rhoDot [m/s]")
    % nexttile
    %     hold on; grid on;
    %     for kk = 1:length(stations)
    %         plot(tMeas, rad2deg(stations(kk).elAngle), 'Color', colors(kk), 'Marker', 'square', 'MarkerSize', 2);
    %         plot(t_nom(k), rad2deg(stations(kk).elAngle(k)), 'Color', colors(kk), 'Marker', '.', 'MarkerSize', 15);
    %     end
    %     yline(rad2deg(elMask),'k--')
    %     xlabel("Time [sec]"); ylabel("Elevation Angle [deg]")
    legend(ax, labels, 'location', 'bestoutside')

        % Update figure and make movie vector
    drawnow
    movieVector = [movieVector; getframe(fig)];

end

    %% Make movie
if saveMovie
    movie = VideoWriter(movieTitle,'MPEG-4');
    movie.FrameRate = 15;
    
        % Open the VideoWriter object, write the movie, and close the file
    open(movie);
    writeVideo(movie, movieVector);
    close(movie);
end
