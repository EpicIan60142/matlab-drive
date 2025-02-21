function fig = plotStateError(t_state, stateError, t_sigma, sigma, boundLevel, titleText, xLabel, yLabel)
% Function that plots state errors from a filter estimate with specified
% sigma bounds
%   Inputs:
%       - t_state: Time vector corresponding to the state errors
%       - stateError: State error matrix, organized as a stack of state
%                     error row vectors - one row per entry in t and one
%                     column per state
%       - t_sigma: Time vector corresponding to the state error uncertainty
%       - sigma: Error uncertainty corresponding to the state errors,
%                organized as a stack of row vectors - one row per entry in
%                t and one column per state
%       - boundLevel: How many multiples of sigma to plot around the state
%                     errors
%       - titleText: Title to be displayed at the top of the plot as a 
%                    string
%       - xLabel: Label of the x-axis as a string
%       - yLabel: String vector of y-axis labels, one for each state in
%                 stateError
%                 Ex. yLabel = ["X Error [km]", "Y Error [km/s]", ...]
%   Outputs:
%       - fig: Graphical handle for the state error figure
%   
%   By: Ian Faber, 02/01/2025
%

fig = figure; tl = tiledlayout(2,3); ax = [];
title(tl, titleText)
for k = 1:size(stateError,2)
    nt = nexttile;
        ax = [ax; nt];
        hold on; grid on;
        plot(t_state, stateError(:,k),'.');
        if ~isempty(sigma)
            plot(t_sigma, -boundLevel*sigma(:,k), 'k--')
            plot(t_sigma, boundLevel*sigma(:,k), 'k--')
        end
        xlabel(xLabel); ylabel(yLabel(k))
        % ylim([-4*max(sigma(:,k)) 4*max(sigma(:,k))]);
end
linkaxes(ax, 'x');

end