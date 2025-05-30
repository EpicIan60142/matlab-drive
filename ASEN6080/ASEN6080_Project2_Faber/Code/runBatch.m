function batchRun = runBatch(X0, x0, P0, pConst, scConst, stations, X_ref, tspan, dt, opt, numIter, measInclude, plot, kEnd)
% Function that runs a batch filter on the given data for a stat OD
% problem. Iterates batch runs until convergence, or the maximum number of
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
%       - stations: Stations structure as defined by makeSations.m
%       - X_ref: Reference orbit from truth data
%       - tspan: Time span of the reference orbit
%       - dt: Output timestep corresponding to tspan
%       - opt: ode45 settings used to generate X_batchFilt from x0Est -
%              should match the settings used to make X_ref
%       - numIter: Max number of times to iterate the filter
%       - measInclude: Boolean array that indicates which measurements to
%                      include in the state estimate. If empty, defaults to
%                      true(1,2) to include both range and range rate
%   - Outputs:
%       - batchRun: Batch filter results structure organized as follows:
%           - batchOut: Filter output structure as defined in BatchFilter.m
%           - t_batchFilt: Propagated batch filter time vector
%           - X_batchFilt: Propagated batch filter state estimates,
%                          organized as follows:
%                          [ 
%                               [X, Y, Z, Xdot, Ydot, Zdot]; 
%                               [X, Y, Z, Xdot, Ydot, Zdot]; ...
%                          ]
%           - P_batch: Propagated batch filter covariance matrices
%           - RMS_prefit_batch: Prefit residual RMS errors for each run 
%                               of the batch filter, one run's error per 
%                               row
%           - RMS_postfit_batch: Postfit residual RMS errors for each run 
%                                of the batch filter, one run's error per 
%                                row
%           - RMS_state_comp_batch: Component-wise state RMS errors after 
%                                   the filter has converged, one row per 
%                                   component
%           - RMS_state_full_batch: Full state RMS error after the filter
%                                   has converged
%           - fig_BatchPreRes: Batch prefit residual plot figure handle
%           - fig_BatchPostRes: Batch postfit residual plot figure handle
%           - fig_BatchPTrace: Batch covariance trace plot figure handle
%           - fig_CovEllipsoids: Covariance ellipsoids plot figure handles
%           - fig_BatchError: Batch filter state error plot figure handle
%
%   By: Ian Faber, 02/03/2025
%

    %% Initialize batch filter
X0_batch = X0;
x0_batch = x0;
P0_batch = P0;

if isempty(measInclude)
    measInclude = true(1,2);
end

if any(~plot)
    fig_BatchPreRes = [];
    fig_BatchPostRes = [];
    fig_BatchPTrace = [];
    fig_CovEllipsoids = [];
    fig_BatchError = [];
end

fprintf("\n\tRunning Batch Filter:\n")

    %% Iterate Batch Filter until residual RMS doesn't change
RMS_prefit_batch = 1e99;
RMS_postfit_batch = 1e99; % Start RMS with a bogus value
batchTolerance = 1e-3; % any ratio less than this is considered converged
maxBatchRuns = numIter; % Cap number of runs
batchRuns = 0;
fig_BatchPreRes = [];
fig_BatchPostRes = [];
k = 2; % Start counter with bogus value
while batchRuns < maxBatchRuns
        % Run batch
    batchOut = BatchFilter(X0_batch, stations, pConst, scConst, P0_batch, x0_batch, dt, kEnd);

        % Extract batch data
    x0Est_batch = batchOut.x0Est;
    P0Est_batch = batchOut.P0Est;
    prefit_res_batch = batchOut.prefit_res;
    postfit_res_batch = batchOut.postfit_res;
    t_batch = batchOut.t;
    statVis_batch = batchOut.statVis;
    
         % Find residual RMS errors
    rms = calcResidualRMS(prefit_res_batch, stations, statVis_batch, measInclude);
    RMS_prefit_batch = [RMS_prefit_batch; rms];

    rms = calcResidualRMS(postfit_res_batch, stations, statVis_batch, measInclude);
    RMS_postfit_batch = [RMS_postfit_batch; rms];
    
        % Plot residuals
    % titleText = sprintf("Batch Filter Pre-Fit Residuals - Run %.0f", k-1); 
    % xLabel = "Time [sec]"; 
    % yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
    % colors = ['b', 'r'];
    
    % fig_BatchPreRes = [fig_BatchPreRes; plotResiduals(t_batch, prefit_res_batch, titleText, xLabel, yLabel, colors)];

    % titleText = sprintf("Batch Filter Post-Fit Residuals - Run %.0f", k-1); 
    % xLabel = "Time [sec]"; 
    % yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
    % colors = ['b', 'r'];
    % 
    % fig_BatchPostRes = [fig_BatchPostRes; plotResiduals(t_batch, postfit_res_batch, titleText, xLabel, yLabel, colors)];

        % Determine if another run is needed via percent change
    if (abs((RMS_postfit_batch(k) - RMS_postfit_batch(k-1))/RMS_postfit_batch(k-1)) > batchTolerance)
            % Update initial state and perturbations for next run
        X0_batch = X0_batch + x0Est_batch;
        x0_batch = x0_batch - x0Est_batch;

            % Update counters
        k = k + 1;
        batchRuns = batchRuns + 1;

            % Write to console
        if batchRuns < maxBatchRuns
            fprintf("Prefit RMS: %.4f, Postfit RMS: %.4f. Iterating Batch. Runs so far: %.0f\n", RMS_prefit_batch(k-1), RMS_postfit_batch(k-1), batchRuns)
        else
            fprintf("Prefit RMS: %.4f, Postfit RMS: %.4f. Hit max Batch iterations. Runs so far: %.0f\n", RMS_prefit_batch(k-1), RMS_postfit_batch(k-1), batchRuns)
        end
    else
        break;
    end

end
RMS_prefit_batch = RMS_prefit_batch(2:end); % Get rid of bogus starting RMS value
RMS_postfit_batch = RMS_postfit_batch(2:end); % Get rid of bogus starting RMS value

if batchRuns < maxBatchRuns
    fprintf("Final prefit RMS: %.4f. Converged after %.0f runs\n", RMS_prefit_batch(end), batchRuns)
    fprintf("Final postfit RMS: %.4f. Converged after %.0f runs\n", RMS_postfit_batch(end), batchRuns)
else
    fprintf("Final prefit RMS: %.4f. Hit maximum number of %.0f runs\n", RMS_prefit_batch(end), maxBatchRuns)
    fprintf("Final postfit RMS: %.4f. Hit maximum number of %.0f runs\n", RMS_postfit_batch(end), maxBatchRuns)
end

    %% Plot residuals and covariance trace
if plot(1)
            % Residuals
    titleText = sprintf("Batch Filter Pre-Fit Residuals - Run %.0f", batchRuns); 
    xLabel = "Time [sec]";
    yLabel = [];
    colors = [];
    if measInclude(1)
        yLabel = [yLabel; "Range Residuals [m]"];
        colors = [colors; 'b'];
    end
    if measInclude(2)
        yLabel = [yLabel;"Range-Rate Residuals [m/s]"];
        colors = [colors; 'r'];
    end
    % yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
    % colors = ['b', 'r'];
    fig_BatchPreRes = plotResiduals(t_batch, prefit_res_batch, titleText, xLabel, yLabel, colors);
end
if plot(2)
    titleText = sprintf("Batch Filter Post-Fit Residuals - Run %.0f", batchRuns); 
    xLabel = "Time [sec]"; 
    yLabel = [];
    colors = [];
    if measInclude(1)
        yLabel = [yLabel; "Range Residuals [m]"];
        colors = [colors; 'b'];
    end
    if measInclude(2)
        yLabel = [yLabel;"Range-Rate Residuals [m/s]"];
        colors = [colors; 'r'];
    end
    % yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
    % colors = ['b', 'r'];
    fig_BatchPostRes = plotResiduals(t_batch, postfit_res_batch, titleText, xLabel, yLabel, colors);
end

Phi = batchOut.Phi;
P_batch = [];
P_plot = [];
% t_P = [];
for k = 1:length(Phi)
    P = Phi{k}*P0Est_batch*Phi{k}';
    P_batch = [P_batch; {P}];
    if any(Phi{k,2} == t_batch) 
        P_plot = [P_plot; {P}];
    end
    % t_P = [t_P; Phi{k,2}];
end

        % Covariance trace
if plot(3)
    titleText = sprintf("Batch Filter R Covariance Trace - Run %.0f", batchRuns); 
    xLabel = "Time [sec]"; 
    yLabel = "trace(P)";
    colors = ['b', 'r'];
    elements = 1:3; % Only plot trace of position
    fig_BatchPTrace = plotPTrace(t_batch, P_plot, elements, titleText, xLabel, yLabel, colors);
end
if plot(4)
    titleText = sprintf("Batch Filter V Covariance Trace - Run %.0f", batchRuns); 
    xLabel = "Time [sec]"; 
    yLabel = "trace(P)";
    colors = ['b', 'r'];
    elements = 4:6; % Only plot trace of velocity
    fig_BatchPTrace = [fig_BatchPTrace; plotPTrace(t_batch, P_plot, elements, titleText, xLabel, yLabel, colors)];
end
    %% Repropagate orbit with the new X0
[t_batchFilt, X_batchFilt] = ode45(@(t,X)orbitEOM_MuSunSRP(t,X,pConst, scConst), tspan, X0_batch, opt);

    %% Calculate state error and uncertainty
stateError_batch = X_batchFilt - X_ref;

sigma_batch = [];
t_sigma = [];
for k = 1:size(batchOut.Phi,1)
    sigPart = [];
    Phi = batchOut.Phi{k,1};
    t_sig = batchOut.Phi{k,2};
    P = Phi*P0Est_batch*Phi';
    for kk = 1:size(P, 1)
        sigPart = [sigPart; sqrt(P(kk,kk))];
    end
    t_sigma = [t_sigma; t_sig];
    sigma_batch = [sigma_batch; sigPart'];
end

    %% Find state RMS error: component-wise and state-wise
[RMS_state_comp_batch, RMS_state_full_batch] = calcStateErrorRMS(stateError_batch);

    %% Plot covariance ellipsoids
fig_CovEllipsoids = [];
P_end = P_batch{end};

if plot(5)
    elements = 1:3;
    P_pos = P_end(elements, elements); mu = X_batchFilt(end,elements)';
    titleText = sprintf("Final Batch Filter Position Covariance Ellipsoid, t = %.3f sec\n\\mu = [%.3e, %.3e, %.3e]^T km\n\\sigma_X = %.3e km, \\sigma_Y = %.3e km, \\sigma_Z = %.3e km", t_batchFilt(end), mu, sqrt(P_pos(1,1)), sqrt(P_pos(2,2)), sqrt(P_pos(3,3)));
    xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; cBarText = "||R - \mu||_2 [km]";
    fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_pos, mu, titleText, xLabel, yLabel, zLabel, cBarText)];
end

if plot(6)
    elements = 4:6;
    P_vel = P_end(elements, elements); mu = X_batchFilt(end,elements)';
    titleText = sprintf("Final Batch Filter Velocity Covariance Ellipsoid, t = %.3f sec\n\\mu = [%.3e, %.3e, %.3e]^T km/s\n\\sigma_{Xdot} = %.3e km/s, \\sigma_{Ydot} = %.3e km/s, \\sigma_{Zdot} = %.3e km/s", t_batchFilt(end), mu, sqrt(P_vel(1,1)), sqrt(P_vel(2,2)), sqrt(P_vel(3,3)));
    xLabel = "Xdot [km/s]"; yLabel = "Ydot [km/s]"; zLabel = "Zdot [km/s]"; cBarText = "||V - \mu||_2 [km/s]";
    fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_vel, mu, titleText, xLabel, yLabel, zLabel, cBarText)];
end

    %% Create state error plots
if plot(7) || plot(8)
    boundLevel = 3; % Plot +/- boundLevel*sigma around state errors
    titleText = "Batch Filter Estimated State Error (X_{filt} - X_{ref})";
    xLabel = "Time [sec]";
    yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
              "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]", ...
              "C_R error [n.d.]"];
    
    fig_BatchError = plotStateError(t_batchFilt, stateError_batch, t_sigma, sigma_batch, boundLevel, titleText, xLabel, yLabel);
end

    %% Assign output
batchRun = struct("batchOut", batchOut, "t_batchFilt", t_batchFilt, "X_batchFilt", X_batchFilt, ...
                  "P_batch", {P_batch}, "RMS_prefit_batch", RMS_prefit_batch, ...
                  "RMS_postfit_batch", RMS_postfit_batch, "RMS_state_comp_batch", RMS_state_comp_batch, ...
                  "RMS_state_full_batch", RMS_state_full_batch, "fig_BatchPreRes", fig_BatchPreRes, ...
                  "fig_BatchPostRes", fig_BatchPostRes, "fig_BatchPTrace", fig_BatchPTrace, ...
                  "fig_CovEllipsoids", fig_CovEllipsoids, "fig_BatchError", fig_BatchError ...
                  );

end