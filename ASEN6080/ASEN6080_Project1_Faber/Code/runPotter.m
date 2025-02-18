function PotterRun = runPotter(X0, x0, P0, pConst, scConst, stations, numIter)
% Function that runs a Potter algorithm on the given data for a stat OD 
% problem. Iterates Potter runs until convergence, or the maximum number of
% iterations is met.
%   - Inputs: 
%       - X0: Initial cartesian state from orbital elements, organized as 
%             [X0; Y0; Z0; Xdot0; Ydot0; Zdot0]
%       - x0: Initial state deviation estimate for the filter to use,
%             organized as
%             [x0; y0; z0; xDot0; yDot0; zDot0]
%       - P0: Initial state covariance matrix for the filter to use
%       - pConst: Planetary constants structure as defined by
%                 getPlanetConst.m
%       - scConst: Spacecraft constants structure as defined by
%                  getSCConst.m
%       - stations: Stations structure as defined by makeSations.m
%       - X_ref: Reference orbit from truth data
%       - t_ref: Reference time vector from truth data\
%       - numIter: Max number of times to iterate the filter
%   - Outputs:
%       - PotterRun: Potter results structure organized as follows:
%           - PotterOut: Filter output structure as defined in Potter.m
%           - t_Potter: Extracted Potter time vector
%           - X_Potter: Extracted Potter state estimates, organized as 
%                       follows:
%                       [ 
%                           [X; Y; Z; Xdot; Ydot; Zdot], 
%                           [X; Y; Z; Xdot; Ydot; Zdot], ...
%                       ]
%           - RMS_prefit_Potter: Prefit residual RMS errors for each run 
%                                of the Potter, one run's error per row
%           - RMS_postfit_Potter: Postfit residual RMS errors for each run 
%                                 of the Potter, one run's error per row
%           - RMS_state_comp_Potter: Component-wise state RMS errors after 
%                                    the filter has converged, one row per 
%                                    component
%           - RMS_state_full_Potter: Full state RMS error after the filter
%                                    has converged
%           - fig_PotterPreRes: Potter prefit residual plot figure handle
%           - fig_PotterPostRes: Potter postfit residual plot figure handle
%           - fig_PotterPTrace: Potter final P trace plot figure handle
%           - fig_CovEllipsoids: Covariance Ellipsoids plot figure handles
%           - fig_PotterError: Potter state error plot figure handle
%
%   By: Ian Faber, 02/03/2025
%
    
    %% Initialize Potter
X0_Potter = X0;
x0_Potter = x0;
P0_Potter = P0;
fprintf("\n\tRunning Potter:\n")

    %% Iterate Potter until residual RMS doesn't change
RMS_prefit_Potter = 1e99;
RMS_postfit_Potter = 1e99; % Start RMS with a bogus value
PotterTolerance = 1e-3; % any ratio less than this is considered converged
maxPotterRuns = numIter; % Cap number of runs
PotterRuns = 0;
% fig_PotterPreRes = [];
% fig_PotterPostRes = [];
k = 2; % Start counter with bogus value
while PotterRuns < maxPotterRuns
        % Run Potter
    PotterOut = Potter(X0_Potter, stations, pConst, scConst, P0_Potter, x0_Potter);

        % Extract Potter data
    xEst_Potter = PotterOut.xEst;
    PEst_Potter = PotterOut.PEst;
    prefit_res_Potter = PotterOut.prefit_res;
    postfit_res_Potter = PotterOut.postfit_res;
    t_Potter = PotterOut.t;
    statVis_Potter = PotterOut.statVis;
    X_Potter = PotterOut.XEst;
    Phi_full_Potter = PotterOut.Phi{end};
    
        % Calculate residual RMS errors
    rms = calcResidualRMS(prefit_res_Potter, stations, statVis_Potter, true(1,2));
    RMS_prefit_Potter = [RMS_prefit_Potter; rms];
    
    rms = calcResidualRMS(postfit_res_Potter, stations, statVis_Potter, true(1,2));
    RMS_postfit_Potter = [RMS_postfit_Potter; rms];

    %     % Plot residuals
    % titleText = sprintf("Potter Pre-Fit Residuals - Run %.0f", k-1); 
    % xLabel = "Time [sec]"; 
    % yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
    % colors = ['b', 'r'];
    % fig_PotterPreRes = [fig_PotterPreRes; plotResiduals(t_Potter, prefit_res_Potter, titleText, xLabel, yLabel, colors)];
    % 
    % titleText = sprintf("Potter Post-Fit Residuals - Run %.0f", k-1); 
    % xLabel = "Time [sec]"; 
    % yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
    % colors = ['b', 'r'];
    % fig_PotterPostRes = [fig_PotterPostRes; plotResiduals(t_Potter, postfit_res_Potter, titleText, xLabel, yLabel, colors)];

        % Determine if another run is needed via percent change
    if (abs((RMS_postfit_Potter(k) - RMS_postfit_Potter(k-1))/RMS_postfit_Potter(k-1)) > PotterTolerance)
            % Update initial state and perturbations for next run
        x0_Potter = (Phi_full_Potter^-1)*xEst_Potter(:,end);
        X0_Potter = X0_Potter + x0_Potter;
        
            % Update counters
        k = k + 1;
        PotterRuns = PotterRuns + 1;

            % Write to console
        if PotterRuns < maxPotterRuns
            fprintf("Prefit RMS: %.4f, Postfit RMS: %.4f. Iterating Potter. Runs so far: %.0f\n", RMS_prefit_Potter(k-1), RMS_postfit_Potter(k-1), PotterRuns)
        else
            fprintf("Prefit RMS: %.4f, Postfit RMS: %.4f. Hit max Potter iterations. Runs so far: %.0f\n", RMS_prefit_Potter(k-1), RMS_postfit_Potter(k-1), PotterRuns)
        end
    else
        break;
    end
end
RMS_prefit_Potter = RMS_prefit_Potter(2:end); % Get rid of bogus starting RMS value
RMS_postfit_Potter = RMS_postfit_Potter(2:end); % Get rid of bogus starting RMS value

if PotterRuns < maxPotterRuns
    fprintf("Final prefit RMS: %.4f. Converged after %.0f runs\n", RMS_prefit_Potter(end), PotterRuns)
    fprintf("Final postfit RMS: %.4f. Converged after %.0f runs\n", RMS_postfit_Potter(end), PotterRuns)
else
    fprintf("Final prefit RMS: %.4f. Hit maximum number of %.0f runs\n", RMS_prefit_Potter(end), PotterRuns)
    fprintf("Final postfit RMS: %.4f. Hit maximum number of %.0f runs\n", RMS_postfit_Potter(end), maxPotterRuns)
end

    %% Plot residuals and covariance trace
        % Residuals
titleText = sprintf("Potter Pre-Fit Residuals - Run %.0f", PotterRuns); 
xLabel = "Time [sec]"; 
yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
colors = ['b', 'r'];
fig_PotterPreRes = plotResiduals(t_Potter, prefit_res_Potter, titleText, xLabel, yLabel, colors);

titleText = sprintf("Potter Post-Fit Residuals - Run %.0f", PotterRuns); 
xLabel = "Time [sec]"; 
yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
colors = ['b', 'r'];
fig_PotterPostRes = plotResiduals(t_Potter, postfit_res_Potter, titleText, xLabel, yLabel, colors);
        
        % Covariance trace
titleText = sprintf("Potter R Covariance Trace - Run %.0f", PotterRuns); 
xLabel = "Time [sec]"; 
yLabel = "trace(P)";
colors = ['b', 'r'];
elements = 1:3; % Only plot trace of position
fig_PotterPTrace = plotPTrace(t_Potter, PotterOut.PEst, elements, titleText, xLabel, yLabel, colors);

titleText = sprintf("Potter V Covariance Trace - Run %.0f", PotterRuns); 
xLabel = "Time [sec]"; 
yLabel = "trace(P)";
colors = ['b', 'r'];
elements = 4:6; % Only plot trace of velocity
fig_PotterPTrace = [fig_PotterPTrace; plotPTrace(t_Potter, PotterOut.PEst, elements, titleText, xLabel, yLabel, colors)];

    %% Plot covariance ellipsoids
fig_CovEllipsoids = [];
P_end = PEst_Potter{end};

elements = 1:3;
P_pos = P_end(elements, elements); mu = X_Potter(elements,end);
titleText = sprintf("Final Potter Position Covariance Ellipsoid\n\\mu = [%.3e, %.3e, %.3e]^T m\n\\sigma_X = %.3e m, \\sigma_Y = %.3e m, \\sigma_Z = %.3e m", mu, sqrt(P_pos(1,1)), sqrt(P_pos(2,2)), sqrt(P_pos(3,3)));
xLabel = "X [m]"; yLabel = "Y [m]"; zLabel = "Z [m]"; cBarText = "||R - \mu||_2 [m]";
fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_pos, mu, titleText, xLabel, yLabel, zLabel, cBarText)];

elements = 4:6;
P_vel = P_end(elements, elements); mu = X_Potter(elements,end);
titleText = sprintf("Final Potter Velocity Covariance Ellipsoid\n\\mu = [%.3e, %.3e, %.3e]^T m/s\n\\sigma_{Xdot} = %.3e m/s, \\sigma_{Ydot} = %.3e m/s, \\sigma_{Zdot} = %.3e m/s", mu, sqrt(P_vel(1,1)), sqrt(P_vel(2,2)), sqrt(P_vel(3,3)));
xLabel = "Xdot [m/s]"; yLabel = "Ydot [m/s]"; zLabel = "Zdot [m/s]"; cBarText = "||V - \mu||_2 [m/s]";
fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_vel, mu, titleText, xLabel, yLabel, zLabel, cBarText)];

    %% Calculate relative state and uncertainty
Phi = PotterOut.Phi;
relState_Potter = [];
for k = 1:length(Phi)
    dX = Phi{k}*x0_Potter - PotterOut.xEst(:,k);
    relState_Potter = [relState_Potter, dX];
end

sigma_Potter = [];
for k = 1:length(PEst_Potter)
    P = PEst_Potter{k};

    sigPart = [];
    for kk = 1:size(P,1)
        sigPart = [sigPart, sqrt(P(kk,kk))];
    end

    sigma_Potter = [sigma_Potter; sigPart];
end

    %% Find state RMS error: component-wise and state-wise
[RMS_state_comp_Potter, RMS_state_full_Potter] = calcStateErrorRMS(relState_Potter);

    %% Create state error plots
boundLevel = 3; % Plot +/- boundLevel*sigma around state errors
titleText = "Potter Relative State (\Deltax_{Potter} = \phi(t,xHat_0,t_0) - xHat(t))";
xLabel = "Time [sec]";
yLabel = ["X error [m]", "Y error [m]", "Z error [m]", ...
          "Xdot error [m/s]", "Ydot error [m/s]", "Zdot error [m/s]", ...
          "\mu error [m^3/s^2]", "J_2 error", "Cd error", ...
          "X_{s,1} error [m]", "Y_{s,1} error [m]", "Z_{s,1} error [m]", ...
          "X_{s,2} error [m]", "Y_{s,2} error [m]", "Z_{s,2} error [m]", ...
          "X_{s,3} error [m]", "Y_{s,3} error [m]", "Z_{s,3} error [m]"];

fig_PotterError = plotStateError(t_Potter, relState_Potter', [], [], boundLevel, titleText, xLabel, yLabel);

    %% Assign output
PotterRun = struct("PotterOut", PotterOut, "t_Potter", t_Potter, "X_Potter", X_Potter, ...
                   "RMS_prefit_Potter", RMS_prefit_Potter, "RMS_postfit_Potter", RMS_postfit_Potter, ...
                   "RMS_state_comp_Potter", RMS_state_comp_Potter, "RMS_state_full_Potter", RMS_state_full_Potter, ...
                   "fig_PotterPreRes", fig_PotterPreRes, "fig_PotterPostRes", fig_PotterPostRes, ...
                   "fig_PotterPTrace", fig_PotterPTrace, "fig_CovEllipsoids", fig_CovEllipsoids, ...
                   "fig_PotterError", fig_PotterError);

end