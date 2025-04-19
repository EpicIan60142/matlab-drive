function fig = cornerPlot(propTraj, monteTraj, P, titleText)
% Function that makes a corner plot of a set of monte carlo runs at a given
% time. 
%   Inputs:
%       -propTraj: Propagated trajectory at the specified time, i.e. via 
%                  nonlinear/LKF/UKF propagation, organized as follows:
%                  [X; Y; Z; Xdot; Ydot; Zdot]
%       - monteTraj: Set of N monte carlo trajectories at the specified
%                    time, given a random initial deviation from nomTraj,
%                    organized as follows:
%                    [monteTraj_1, monteTraj_2, ..., monteTraj_N], where
%                    monteTraj_N = [X; Y; Z; Xdot; Ydot; Zdot]
%       - titleText: Title text for this corner plot
%   Outputs:
%       - fig: Corner plot figure handle
%
%   By: Ian Faber, 04/15/2025   
%

    % Calculate mean of monte carlo trajectories
mu = mean(monteTraj,2);

% if isempty(sigma)
%     sigma = std(monteTraj,0,2);
% end
    % Construct covariance matrix from monte carlo trajectories
if isempty(P)
    P = cov(monteTraj'); % Need input to be structured with observations (runs) as rows and variables as columns
end
    % Extract valid covariances for corner plotting
Ps = {}; % Cell array for holding valid covariances
stateIdx = []; % Valid states for the saved covariances, i.e. whether they belong to X, ZYdot, etc.
for vert = 1:6 % Vertical subplot coordinate X -> Zdot = 1 -> 6, move top -> down
    for horiz = 1:6 % Horizontal subplot coordinate X -> Zdot = 1 -> 6, move left -> right
        if horiz > vert % Only interested in lower triangular covariances, i.e. horizontal coordinate smaller than vertical coordinate
            continue;
        else
            if horiz == vert
                Pcorner = P(horiz, horiz);
                stateIdx = [stateIdx, [horiz; vert]];
            else
                Pcorner = P([horiz,vert],[horiz,vert]);
                stateIdx = [stateIdx, [horiz; vert]];
            end
            Ps = [Ps; {Pcorner}];
        end
    end
end

    % Create helper variables for ellipses
theta = linspace(0, 2*pi, 100);
circle.x = cos(theta);
circle.y = sin(theta);

    % Create corner plots
validTiles = [1, 7, 8, 13, 14, 15, 19, 20, 21, 22, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36]; % Valid plotting tiles
labels = ["X [km]", "Y [km]", "Z [km]", "Xdot [km/s]", "Ydot [km/s]", "Zdot [km/s]"];
fig = figure; tl = tiledlayout(6,6); tiles = [];
title(tl, titleText)
for k = 1:length(validTiles)
    nt = nexttile(validTiles(k));
        hold on; grid on;
        if size(Ps{k},1) == 1 && size(Ps{k},2) == 1 % This is a single variable variance - need histogram
                % Create pdf
            % x = linspace(min(monteTraj(stateIdx(1,k),:)),max(monteTraj(stateIdx(1,k),:)),100);
            sig = sqrt(Ps{k});
            x = linspace(mu(stateIdx(1,k)) - 3*sig, mu(stateIdx(1,k)) + 3*sig, 1000);
            px = normpdf(x, mu(stateIdx(1,k)), sqrt(Ps{k}));%sigma(stateIdx(1,k)));
                % Plot histogram and pdf
            histogram(monteTraj(stateIdx(1,k),:),'FaceColor','b','Normalization','pdf');
            pdf = plot(x, px, 'm--', 'LineWidth', 2);
                % Plot propagated state coordinate
            nomSingle = xline(propTraj(stateIdx(1,k)), 'k--', 'LineWidth', 2);
                % Labels
            xlabel(labels(stateIdx(1,k)));
        else % This is a covariance - need ellipse
                % Find eigenvectors and values
            [eigVec, Lambda] = eig(Ps{k});
            Lambda = chol(Lambda);
                % Find ellipse axis scales
            scales = [];
            for kk = 1:size(Lambda,1)
                scales = [scales; Lambda(kk,kk)];
            end
                % Create 1, 2, and 3 sigma ellipses
            ellipse_1sig = eigVec*[scales(1)*circle.x; scales(2)*circle.y] + mu([stateIdx(1,k), stateIdx(2,k)]);
            ellipse_2sig = eigVec*[2*scales(1)*circle.x; 2*scales(2)*circle.y] + mu([stateIdx(1,k), stateIdx(2,k)]);
            ellipse_3sig = eigVec*[3*scales(1)*circle.x; 3*scales(2)*circle.y] + mu([stateIdx(1,k), stateIdx(2,k)]);
                % Plot ellipses
            Sig_1 = plot(ellipse_1sig(1,:), ellipse_1sig(2,:), 'b--');
            Sig_2 = plot(ellipse_2sig(1,:), ellipse_2sig(2,:), 'r--');
            Sig_3 = plot(ellipse_3sig(1,:), ellipse_3sig(2,:), 'g--');
                % Plot all monte carlo points
            montePoint = scatter(monteTraj(stateIdx(1,k),:), monteTraj(stateIdx(2,k),:), 3, 'black', 'filled', 'o', 'MarkerFaceAlpha', 0.25);
                % Plot mean coordinate
            meanPoint = scatter(mu(stateIdx(1,k)), mu(stateIdx(2,k)), 25, 'black', 'filled' , 'square');
                % Plot propagated coordinate
            nomPoint = scatter(propTraj(stateIdx(1,k)), propTraj(stateIdx(2,k)), 25, 'black', 'filled', '^');
                % Labels
            xlabel(labels(stateIdx(1,k))); ylabel(labels(stateIdx(2,k)));
        end
    tiles = [tiles; nt];
end

lgnd = legend([pdf, nomSingle, Sig_1, Sig_2, Sig_3, montePoint, meanPoint, nomPoint], ...
              ["State Component pdf", "State Component Nominal value", "1\sigma Ellipse", "2\sigma Ellipse", "3\sigma Ellipse", "MC Points", "MC Mean", "Propagated Trajectory Point"], ...
               'Location', 'layout');
lgnd.Layout.Tile = 11;

    % Link x and y axes where appropriate
        % x axes
for k = 1:6
    idx = find(stateIdx(1,:) == k);
    ax = tiles(idx);
    linkaxes(ax, 'x');
end
        % y axes
for k = 1:6
    idy = find(stateIdx(2,:) == k);
    if length(idy) > 1
        ay = tiles(idy(1:end-1)); % don't include histograms, which are the last tiles plotted in each vertical coordinate
    else % only X histogram extracted
        ay = [];
    end

    if ~isempty(ay)
        linkaxes(ay, 'y');
    end

end


end