function LKFRun = runLKF(X0, x0, P0, pConst, stations, X_ref, t_ref, numIter)
% Function that runs an LKF on the given data for a stat OD problem. 
% Iterates LKF runs until convergence, or the maximum number of iterations 
% is met.
%   - Inputs: 
%       - X0: Initial cartesian state from orbital elements, organized as 
%             [X0; Y0; Z0; Xdot0; Ydot0; Zdot0]
%       - x0: Initial state deviation estimate for the filter to use,
%             organized as
%             [x0; y0; z0; xDot0; yDot0; zDot0]
%       - P0: Initial state covariance matrix for the filter to use
%       - pConst: Planetary constants structure as defined by
%                 getPlanetConst.m
%       - stations: Stations structure as defined by makeSations.m
%       - X_ref: Reference orbit from truth data
%       - t_ref: Reference time vector from truth data
%       - numIter: Max number of times to iterate the filter
%   - Outputs:
%       - LKFRun: LKF results structure organized as follows:
%           - LKFOut: Filter output structure as defined in LKF.m
%           - t_LKF: Extracted LKF time vector
%           - X_LKF: Extracted LKF state estimates, organized as follows:
%                    [ 
%                         [X; Y; Z; Xdot; Ydot; Zdot], 
%                         [X; Y; Z; Xdot; Ydot; Zdot], ...
%                    ]
%           - RMS_prefit_LKF: Prefit residual RMS errors for each run 
%                             of the LKF, one run's error per row
%           - RMS_postfit_LKF: Postfit residual RMS errors for each run 
%                              of the LKF, one run's error per row
%           - RMS_state_comp_LKF: Component-wise state RMS errors after 
%                                 the filter has converged, one row per 
%                                 component
%           - RMS_state_full_LKF: Full state RMS error after the filter
%                                 has converged
%           - fig_LKFPreRes: LKF prefit residual plot figure handle
%           - fig_LKFPostRes: LKF postfit residual plot figure handle
%           - fig_LKFPTrace: LKF final P trace plot figure handle
%           - fig_CovEllipsoids: Covariance Ellipsoids plot figure handles
%           - fig_LKFError: LKF state error plot figure handle
%
%   By: Ian Faber, 02/03/2025
%
    
    %% Initialize LKF
X0_LKF = X0;
x0_LKF = x0;
P0_LKF = P0;
fprintf("\n\tRunning LKF:\n")

    %% Iterate LKF until residual RMS doesn't change
RMS_prefit_LKF = 1e99;
RMS_postfit_LKF = 1e99; % Start RMS with a bogus value
LKFTolerance = 1e-3; % any ratio less than this is considered converged
maxLKFRuns = numIter; % Cap number of runs
LKFRuns = 0;
% fig_LKFPreRes = [];
% fig_LKFPostRes = [];
k = 2; % Start counter with bogus value
while LKFRuns < maxLKFRuns
        % Run LKF
    LKFOut = LKF(X0_LKF, stations, pConst, P0_LKF, x0_LKF);

        % Extract LKF data
    xEst_LKF = LKFOut.xEst;
    PEst_LKF = LKFOut.PEst;
    prefit_res_LKF = LKFOut.prefit_res;
    postfit_res_LKF = LKFOut.postfit_res;
    t_LKF = LKFOut.t;
    statVis_LKF = LKFOut.statVis;
    X_LKF = LKFOut.XEst;
    Phi_full_LKF = LKFOut.Phi{end};
    
        % Calculate residual RMS errors
    rms = calcResidualRMS(prefit_res_LKF, stations, statVis_LKF, true(1,2));
    RMS_prefit_LKF = [RMS_prefit_LKF; rms];
    
    rms = calcResidualRMS(postfit_res_LKF, stations, statVis_LKF, true(1,2));
    RMS_postfit_LKF = [RMS_postfit_LKF; rms];

        % Plot residuals
    % titleText = sprintf("LKF Pre-Fit Residuals - Run %.0f", k-1); 
    % xLabel = "Time [sec]"; 
    % yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
    % colors = ['b', 'r'];
    % fig_LKFPreRes = [fig_LKFPreRes; plotResiduals(t_LKF, prefit_res_LKF, titleText, xLabel, yLabel, colors)];
    % 
    % titleText = sprintf("LKF Post-Fit Residuals - Run %.0f", k-1); 
    % xLabel = "Time [sec]"; 
    % yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
    % colors = ['b', 'r'];
    % fig_LKFPostRes = [fig_LKFPostRes; plotResiduals(t_LKF, postfit_res_LKF, titleText, xLabel, yLabel, colors)];

        % Determine if another run is needed via percent change
    if (abs((RMS_postfit_LKF(k) - RMS_postfit_LKF(k-1))/RMS_postfit_LKF(k-1)) > LKFTolerance)
            % Update initial state and perturbations for next run
        x0_LKF = (Phi_full_LKF^-1)*xEst_LKF(:,end);
        X0_LKF = X0_LKF + x0_LKF;

            % Update counters
        k = k + 1;
        LKFRuns = LKFRuns + 1;

            % Write to console
        if LKFRuns < maxLKFRuns
            fprintf("Prefit RMS: %.4f, Postfit RMS: %.4f. Iterating LKF. Runs so far: %.0f\n", RMS_prefit_LKF(k-1), RMS_postfit_LKF(k-1), LKFRuns)
        else
            fprintf("Prefit RMS: %.4f, Postfit RMS: %.4f. Hit max LKF iterations. Runs so far: %.0f\n", RMS_prefit_LKF(k-1), RMS_postfit_LKF(k-1), LKFRuns)
        end
    else
        break;
    end
end
RMS_prefit_LKF = RMS_prefit_LKF(2:end); % Get rid of bogus starting RMS value
RMS_postfit_LKF = RMS_postfit_LKF(2:end); % Get rid of bogus starting RMS value

if LKFRuns < maxLKFRuns
    fprintf("Final prefit RMS: %.4f. Converged after %.0f runs\n", RMS_prefit_LKF(end), LKFRuns)
    fprintf("Final postfit RMS: %.4f. Converged after %.0f runs\n", RMS_postfit_LKF(end), LKFRuns)
else
    fprintf("Final prefit RMS: %.4f. Hit maximum number of %.0f runs\n", RMS_prefit_LKF(end), LKFRuns)
    fprintf("Final postfit RMS: %.4f. Hit maximum number of %.0f runs\n", RMS_postfit_LKF(end), maxLKFRuns)
end

    %% Plot residuals and covariance trace
        % Residuals
titleText = sprintf("LKF Pre-Fit Residuals - Run %.0f", LKFRuns); 
xLabel = "Time [sec]"; 
yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
colors = ['b', 'r'];
fig_LKFPreRes = plotResiduals(t_LKF, prefit_res_LKF, titleText, xLabel, yLabel, colors);

titleText = sprintf("LKF Post-Fit Residuals - Run %.0f", LKFRuns); 
xLabel = "Time [sec]"; 
yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
colors = ['b', 'r'];
fig_LKFPostRes = plotResiduals(t_LKF, postfit_res_LKF, titleText, xLabel, yLabel, colors);
        
        % Covariance trace
titleText = sprintf("LKF R Covariance Trace - Run %.0f", LKFRuns); 
xLabel = "Time [sec]"; 
yLabel = "trace(P)";
colors = ['b', 'r'];
elements = 1:3; % Only plot trace of position
fig_LKFPTrace = plotPTrace(t_LKF, LKFOut.PEst, elements, titleText, xLabel, yLabel, colors);

titleText = sprintf("LKF V Covariance Trace - Run %.0f", LKFRuns); 
xLabel = "Time [sec]"; 
yLabel = "trace(P)";
colors = ['b', 'r'];
elements = 4:6; % Only plot trace of velocity
fig_LKFPTrace = [fig_LKFPTrace; plotPTrace(t_LKF, LKFOut.PEst, elements, titleText, xLabel, yLabel, colors)];

    %% Plot covariance ellipsoids
fig_CovEllipsoids = [];
P_end = PEst_LKF{end};

elements = 1:3;
P_pos = P_end(elements, elements); mu = X_LKF(elements,end);
titleText = sprintf("Final LKF Position Covariance Ellipsoid\n\\mu = [%.3e, %.3e, %.3e]^T m\n\\sigma_X = %.3e m, \\sigma_Y = %.3e m, \\sigma_Z = %.3e m", mu, sqrt(P_pos(1,1)), sqrt(P_pos(2,2)), sqrt(P_pos(3,3)));
xLabel = "X [m]"; yLabel = "Y [m]"; zLabel = "Z [m]"; cBarText = "||R - \mu||_2 [m]";
fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_pos, mu, titleText, xLabel, yLabel, zLabel, cBarText)];

elements = 4:6;
P_vel = P_end(elements, elements); mu = X_LKF(elements,end);
titleText = sprintf("Final LKF Velocity Covariance Ellipsoid\n\\mu = [%.3e, %.3e, %.3e]^T m/s\n\\sigma_{Xdot} = %.3e m/s, \\sigma_{Ydot} = %.3e m/s, \\sigma_{Zdot} = %.3e m/s", mu, sqrt(P_vel(1,1)), sqrt(P_vel(2,2)), sqrt(P_vel(3,3)));
xLabel = "Xdot [m/s]"; yLabel = "Ydot [m/s]"; zLabel = "Zdot [m/s]"; cBarText = "||V - \mu||_2 [m/s]";
fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_vel, mu, titleText, xLabel, yLabel, zLabel, cBarText)];

    %% Extract nominal state from LKF deviation estimates
X_ref_LKF = [];
for k = 1:length(t_LKF)
    X_ref_LKF = [X_ref_LKF; X_ref(t_ref == t_LKF(k),:)];
end

    %% Calculate state error and uncertainty
stateError_LKF = X_LKF' - X_ref_LKF;

sigma_LKF = [];
for k = 1:length(PEst_LKF)
    P = PEst_LKF{k};

    sigPart = [];
    for kk = 1:size(P,1)
        sigPart = [sigPart, sqrt(P(kk,kk))];
    end

    sigma_LKF = [sigma_LKF; sigPart];
end

    %% Find state RMS error: component-wise and state-wise
[RMS_state_comp_LKF, RMS_state_full_LKF] = calcStateErrorRMS(stateError_LKF);

    %% Create state error plots
boundLevel = 3; % Plot +/- boundLevel*sigma around state errors
titleText = "LKF Estimated State Error (X_{filt} - X_{ref})";
xLabel = "Time [sec]";
yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
          "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]"];

fig_LKFError = plotStateError(t_LKF, relState_LKF', [], [], boundLevel, titleText, xLabel, yLabel);
fig_LKFError = plotStateError(t_LKF, relState_LKF', t_LKF, sigma_LKF, boundLevel, titleText, xLabel, yLabel);

    %% Assign output
LKFRun = struct("LKFOut", LKFOut, "t_LKF", t_LKF, "X_LKF", X_LKF, ...
                "RMS_prefit_LKF", RMS_prefit_LKF, "RMS_postfit_LKF", RMS_postfit_LKF, ...
                "RMS_state_comp_LKF", RMS_state_comp_LKF, "RMS_state_full_LKF", RMS_state_full_LKF, ...
                "fig_LKFPreRes", fig_LKFPreRes, "fig_LKFPostRes", fig_LKFPostRes, ...
                "fig_LKFPTrace", fig_LKFPTrace, "fig_CovEllipsoids", fig_CovEllipsoids, ...
                "fig_LKFError", fig_LKFError);

end