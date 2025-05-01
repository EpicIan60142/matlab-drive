function fig_Meas = plotMeasurements(stations, titleText, xLabel, yLabel)
% Function that plots measurements for a given truth data set
%   Inputs: 
%       - stations: stations struct as formatted by makeStations.m
%       - titleText: Title to display on the plot
%       - xLabel: x label for each plot
%       - yLabel: Vector of y label strings for each measurement to plot

fig_Meas = figure; tl = tiledlayout(length(yLabel),1); ax = []; lx = [];
title(tl, titleText)
nt = nexttile;
    hold on; grid on;
    for k = 1:length(stations)
        ax = [ax; plot(stations(k).tMeas, stations(k).rho, '.', 'Color', stations(k).color)];
        label(k) = sprintf("%s",stations(k).id);
    end
    xlabel(xLabel); ylabel(yLabel(1)); lx = [lx; nt];
nt = nexttile;
    hold on; grid on;
    for k = 1:length(stations)
        plot(stations(k).tMeas, stations(k).rhoDot, '.', 'Color', stations(k).color)
    end
    xlabel(xLabel); ylabel(yLabel(2)); lx = [lx; nt];
if length(yLabel) > 2
    nt = nexttile;
        hold on; grid on;
        for k = 1:length(stations)
            plot(stations(k).tMeas, stations(k).elAngle, '.', 'Color', stations(k).color)
        end
        xlabel(xLabel); ylabel(yLabel(3)); lx = [lx; nt];
end
legend(ax, label, 'location', 'bestoutside')
linkaxes(lx, 'x');
drawnow;

end