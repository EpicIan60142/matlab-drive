function fig = plotResiduals(t, residuals, titleText, xLabel, yLabel, colors)
% Function that plots residuals of a filter run
%   Inputs: 
%       - t: Time vector corresponding to the post-fit residuals
%       - residuals: Matrix containing the residuals (pre or post-fit) for
%                    a filter run, where each set of residuals is a row
%                    vector
%                    Ex. residuals = [residualVec_1; residualVec_2]
%       - titleText: Title to be displayed at the top of the plot as a 
%                    string
%       - xLabel: Label of the x-axis as a string
%       - yLabel: String vector of y-axis labels, one for each residual
%                 type in postfit_res
%                 Ex. yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"]
%       - colors: Vector of colors to plot residuals in, specified as color
%                 characters
%                 Ex. colors = ['b', 'r']
%   Outputs:
%       - fig: Graphical handle of the residual plot
%
%   By: Ian Faber, 02/01/2025
%
fig = figure; tl = tiledlayout(size(residuals,1),3); ax = [];
title(tl, titleText)
for k = 1:size(residuals,1)
    nt = nexttile([1 2]); % Span 1 row and 2 columns
        ax = [ax; nt];
        hold on; grid on;
        plot(t, residuals(k,:), '.', 'Color', colors(k))
        xlabel(xLabel); ylabel(yLabel(k))
    hist = nexttile;
        hold on; grid on;
        mu = mean(residuals(k,:)); stdev = std(residuals(k,:));
        infoText = sprintf("\\mu = %.4e \n\\sigma = %.4e", mu, stdev);
        histogram(residuals(k,:), 'FaceColor', colors(k), 'Orientation', 'horizontal')
        yl = ylim; xl = xlim;
        text(0.5*xl(2),0.9*yl(2),infoText, "EdgeColor", 'k', 'BackgroundColor', 'w')
    linkaxes([nt; hist], 'y');
end
linkaxes(ax, 'x')

end