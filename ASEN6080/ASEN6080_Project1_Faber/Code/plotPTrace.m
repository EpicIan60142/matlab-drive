function fig = plotPTrace(t, P, elements, titleText, xLabel, yLabel, colors)
% Function that plots the covariance matrix trace of a filter run
%   Inputs: 
%       - t: Time vector corresponding to the covariance matrices
%       - P: Cell array of covariance matrices, one entry per time in t
%       - elements: Elements of P to include in the trace calculation
%       - titleText: Title to be displayed at the top of the plot as a 
%                    string
%       - xLabel: Label of the x-axis as a string
%       - yLabel: String vector of y-axis labels, one for each residual
%                 type in postfit_res
%                 Ex. yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"]
%       - colors: Vector of colors to plot trace in, specified as color
%                 characters
%                 Ex. colors = ['b', 'r']
%   Outputs:
%       - fig: Graphical handle of the residual plot
%
%   By: Ian Faber, 02/01/2025
%

traces = [];
for k = 1:length(P)
    Pk = P{k};
    traces = [traces; trace(Pk(elements,elements))];
end

fig = figure;
if any(traces<0)
    yyaxis left
end
a = semilogy(t(traces>0), traces(traces>0), '.', 'Color', colors(1));
yLim = ylim;
xlabel(xLabel); ylabel(yLabel);
% hold on; grid on;
if any(traces<0)
    yyaxis right
    negTrace = -traces(traces<0);
    b = semilogy(t(traces<0), negTrace, '.', 'Color', colors(2), 'MarkerSize', 15);
    ylim(yLim)
    xlabel(xLabel), ylabel("-trace(P)")
    legend([a,b], ["trace(P)>0","trace(P)<0"])
end
grid on;
title(titleText)


end