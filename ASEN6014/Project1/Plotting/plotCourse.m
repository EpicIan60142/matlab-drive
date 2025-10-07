function fig = plotCourse(rings, figNum, titleText, xLabel, yLabel, zLabel)
% Function that plots a generated cubesat race course
%   Inputs:
%       - rings: Vector of ring structures for the intermediate rings of
%                the course
%       - figNum: Figure number for the course figure
%       - titleText: String specifying the title for the figure
%       - xLabel: String specifying the x axis for the 3D plot
%       - yLabel: String specifying the y axis for the 3D plot
%       - zLabel: String specifying the z axis for the 3D plot
%   Outputs:
%       - fig: Course figure handle
%
%   By: Ian Faber, 05/05/2025
%

    %% Normal vector scale
normalScale = 0.1;

    %% Make figure
fig = figure(figNum);
fig.WindowState = "maximized";
% tl = tiledlayout(1,2);
% title(tl, titleText);

    %% Plot course
% nexttile(1)
    title(titleText);
    hold on; grid on; axis equal
    
        % Plot intermediate rings
    for k = 2:length(rings)-1
        ring = scatter3(rings(k).center(1), rings(k).center(2), rings(k).center(3), 20, k, 'filled');
        normal = quiver3(rings(k).center(1), rings(k).center(2), rings(k).center(3), rings(k).normal(1), rings(k).normal(2), rings(k).normal(3), normalScale, 'filled', 'k-');
        plotRing(rings(k), 'k-');
    end
        
        % Plot start and end rings
    startRing = rings(1);
    endRing = rings(end);
    cubeStart = plotRing(startRing, 'g-'); cubeStart.LineWidth = 2;
    quiver3(startRing.center(1), startRing.center(2), startRing.center(3), startRing.normal(1), startRing.normal(2), startRing.normal(3), normalScale, 'filled', 'k-')
    cubeEnd = plotRing(endRing, 'r-');
    
        % Plot course origin
    courseCenter = scatter3(0, 0, 0, 20, 'k', 'filled', 'h');
    
        % Labels, colorbar, and view angle
    xlabel(xLabel); ylabel(yLabel); zlabel(zLabel); cBar = colorbar;
    cBar.Label.String = "Ring Number"; cBar.Location = "west";
    colormap("cool"); view([30 35])

    %% Legend
lgnd = legend([cubeStart, cubeEnd, courseCenter, ring, normal], ["CubeSat 3\sigma Starting Ring", "Race Course End Ring", "Race Course Origin", "Course Ring", "Ring Normal Vector"], 'location', 'eastoutside');
% lgnd.Layout.Tile = 2;

end