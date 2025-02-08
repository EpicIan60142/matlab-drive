function fig_Meas = plotMeasurements(stations, titleText, xLabel, yLabel)
% Function that plots measurements for a given truth data set
%   Inputs: 
%       - stations: stations struct as formatted by makeStations.m
%       - titleText: Title to display on the plot
%       - xLabel: x label for each plot
%       - yLabel: Vector of y label strings for each measurement

fig_Meas = figure; tl = tiledlayout(3,1); ax = [];
title(tl, titleText)
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        ax = [ax; plot(stations(k).tMeas, stations(k).rho, '.', 'Color', stations(k).color)];
        label(k) = sprintf("Station %.0f", k);
    end
    xlabel(xLabel); ylabel(yLabel(1));
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        plot(stations(k).tMeas, stations(k).rhoDot, '.', 'Color', stations(k).color)
    end
    xlabel(xLabel); ylabel(yLabel(2));
nexttile
    hold on; grid on;
    for k = 1:length(stations)
        plot(stations(k).tMeas, stations(k).elAngle, '.', 'Color', stations(k).color)
    end
    xlabel(xLabel); ylabel(yLabel(3));
legend(ax, label, 'location', 'bestoutside')

end