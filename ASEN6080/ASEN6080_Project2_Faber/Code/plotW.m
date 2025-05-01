function fig = plotW(t, X, titleText, xLabel, yLabel, colors)
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
    % Pull out w's
W = X(8:10,:);

fig = figure; tl = tiledlayout(3,1); ax = [];
title(tl, titleText)
for k = 1:3
    nt = nexttile;
        ax = [ax; nt];
        hold on; grid on;
        plot(t, W(k,:), '.', 'Color', colors(k))
        xlabel(xLabel); ylabel(yLabel(k))
end
linkaxes(ax, 'x')

end