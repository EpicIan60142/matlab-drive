function fig = plotSegment(cubesat, ring, t, X, u, figNum, titleText, xLabel, yLabel, zLabel)
% Function that plots a single segment of the cubesat race course
%   Inputs:
%       - cubesat: Cubesat structure for the cubesat of interest
%       - ring: The last ring of this segment
%       - t: Time vector for this segment
%       - X: State vector for this segment
%       - figNum: Figure number for this segment
%       - 
%   Outputs:
%       - fig: Figure handle for this segment

    %% Create figure
% if newFigure
fig = figure(figNum); 
fig.WindowState = "maximized";
tl = tiledlayout(5, 7);
title(tl, titleText);
% else
%     fig = figure(figNum);
%     fig.WindowState = "maximized";
%     nexttile(1);
% end

    %% Utilities
marker = cubesat.marker;
markerSize = 10;
color = cubesat.color;
ax = [];

    %% 3 tiles x 3 tiles for 3D view
nexttile([5 4]) 
    hold on; grid on; axis equal
    plotRing(ring.params.lastRing, 'g-');
    plot3(X(1,1), X(1,2), X(1,3), 'k', 'Marker', marker);
    plotRing(ring, 'r-');
    plot3(X(:,1), X(:,2), X(:,3), '--', 'Color', color);
    xlabel(xLabel); ylabel(yLabel); zlabel(zLabel);
    view([30 35]);

    %% Position
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,1), 'Color', color);
    plot(t(1), ring.params.lastRing.center(1), 'g.', 'MarkerSize', markerSize);
    plot(t(end), ring.center(1), 'r.', 'MarkerSize', markerSize);
    xlabel("Time [sec]");
    ylabel("X [m]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,2), 'Color', color);
    plot(t(1), ring.params.lastRing.center(2), 'g.', 'MarkerSize', markerSize);
    plot(t(end), ring.center(2), 'r.', 'MarkerSize', markerSize);
    xlabel("Time [sec]");
    ylabel("Y [m]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,3), 'Color', color);
    plot(t(1), ring.params.lastRing.center(3), 'g.', 'MarkerSize', markerSize);
    plot(t(end), ring.center(3), 'r.', 'MarkerSize', markerSize);
    xlabel("Time [sec]");
    ylabel("Z [m]");

    %% Velocity
vf = vecnorm(X(:,4:6), 2, 2);
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,4), 'Color', color);
    plot(t(1), ring.params.lastRing.normal(1)*vf(1), 'g.', 'MarkerSize', markerSize);
    plot(t(end), ring.normal(1)*vf(end), 'r.', 'MarkerSize', markerSize);
    xlabel("Time [sec]");
    ylabel("Xdot [m/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,5), 'Color', color);
    plot(t(1), ring.params.lastRing.normal(2)*vf(1), 'g.', 'MarkerSize', markerSize);
    plot(t(end), ring.normal(2)*vf(end), 'r.', 'MarkerSize', markerSize);
    xlabel("Time [sec]");
    ylabel("Ydot [m/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,6), 'Color', color);
    plot(t(1), ring.params.lastRing.normal(3)*vf(1), 'g.', 'MarkerSize', markerSize);
    plot(t(end), ring.normal(3)*vf(end), 'r.', 'MarkerSize', markerSize);
    xlabel("Time [sec]");
    ylabel("Zdot [m/s]");

    %% Position Adjoints
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,7), 'Color', color);
    xlabel("Time [sec]");
    ylabel("p_x [m]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,8), 'Color', color);
    xlabel("Time [sec]");
    ylabel("p_y [m]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,9), 'Color', color);
    xlabel("Time [sec]");
    ylabel("p_z [m]");

    %% Velocity Adjoints
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,10), 'Color', color);
    xlabel("Time [sec]");
    ylabel("p_{vx} [m/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,11), 'Color', color);
    xlabel("Time [sec]");
    ylabel("p_{vy} [m/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,12), 'Color', color);
    xlabel("Time [sec]");
    ylabel("p_{vz} [m/s]");

    %% Control
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, u(:,1), 'Color', color);
    xlabel("Time [sec]");
    ylabel("u_x [m/s^2]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, u(:,2), 'Color', color);
    xlabel("Time [sec]");
    ylabel("u_y [m/s^2]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, u(:,3), 'Color', color);
    xlabel("Time [sec]");
    ylabel("u_z [m/s^2]");

    %% Link axes and update figure
linkaxes(ax,'x');
drawnow;
end