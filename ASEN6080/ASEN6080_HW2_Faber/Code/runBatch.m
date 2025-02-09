function batchRun = runBatch(X0, x0, P0, pConst, stations, X_ref, tspan, dt, opt)
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
%       - opt: ode45 settings used to generate X_ref
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
%           - RMS_postfit_batch: Postfit residual RMS errors for each run 
%                                of the batch filter, one run's error per 
%                                row
%           - RMS_state_comp_batch: Component-wise state RMS errors after 
%                                   the filter has converged, one row per 
%                                   component
%           - RMS_state_full_batch: Full state RMS error after the filter
%                                   has converged
%           - fig_BatchPreRes: Array of batch prefit residual plot figure 
%                              handles, one per filter iteration
%           - fig_BatchPostRes: Array of batch postfit residual plot figure 
%                              handles, one per filter iteration
%           - fig_BatchError: Batch filter state error plot figure handle
%
%   By: Ian Faber, 02/03/2025
%

    %% Initialize batch filter
X0_batch = X0;
x0_batch = x0;
P0_batch = P0;

fprintf("\n\tRunning Batch Filter:\n")

    %% Iterate Batch Filter until residual RMS doesn't change
RMS_postfit_batch = 1e99; % Start RMS with a bogus value
batchTolerance = 1e-6; % any ratio less than this is considered converged
maxBatchRuns = 10; % Cap number of runs
batchRuns = 0;
fig_BatchPreRes = [];
fig_BatchPostRes = [];
k = 2; % Start counter with bogus value
while true
        % Run batch
    batchOut = BatchFilter(X0_batch, stations, pConst, P0_batch, x0_batch, dt);

        % Extract batch data
    x0Est_batch = batchOut.x0Est;
    P0Est_batch = batchOut.P0Est;
    prefit_res_batch = batchOut.prefit_res;
    postfit_res_batch = batchOut.postfit_res;
    t_batch = batchOut.t;
    statVis_batch = batchOut.statVis;
    
        % Find residual RMS errors
    rms = calcResidualRMS(postfit_res_batch, stations, statVis_batch);
    RMS_postfit_batch = [RMS_postfit_batch; rms];
    
        % Plot residuals
    titleText = sprintf("Batch Filter Pre-Fit Residuals - Run %.0f", k-1); 
    xLabel = "Time [sec]"; 
    yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
    colors = ['b', 'r'];
    
    fig_BatchPreRes = [fig_BatchPreRes; plotResiduals(t_batch, prefit_res_batch, titleText, xLabel, yLabel, colors)];

    titleText = sprintf("Batch Filter Post-Fit Residuals - Run %.0f", k-1); 
    xLabel = "Time [sec]"; 
    yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
    colors = ['b', 'r'];
    
    fig_BatchPostRes = [fig_BatchPostRes; plotResiduals(t_batch, postfit_res_batch, titleText, xLabel, yLabel, colors)];

        % Determine if another run is needed via percent change
    if (abs((RMS_postfit_batch(k) - RMS_postfit_batch(k-1))/RMS_postfit_batch(k-1)) > batchTolerance) && (batchRuns < maxBatchRuns)
            % Update initial state and perturbations for next run
        X0_batch = X0_batch + x0Est_batch;
        x0_batch = x0_batch - x0Est_batch;

            % Update counters
        k = k + 1;
        batchRuns = k - 1; % Offset current k by 1 to account for bogus starting value

            % Write to console
        fprintf("Postfit RMS: %.4f. Iterating Batch with x0_batch = [%.3e; %.3e; %.3e; %.3e; %.3e; %.3e]. Runs so far: %.0f\n", RMS_postfit_batch(k-1), x0_batch, batchRuns-1)
    else
        break;
    end
end
RMS_postfit_batch = RMS_postfit_batch(2:end); % Get rid of bogus starting RMS value

if batchRuns < maxBatchRuns
    fprintf("Final postfit RMS: %.4f. Converged after %.0f runs\n", RMS_postfit_batch(end), batchRuns)
else
    fprintf("Final postfit RMS: %.4f. Hit maximum number of %.0f runs\n", RMS_postfit_batch(end), maxBatchRuns)
end

    %% Repropagate orbit with the new X0
[t_batchFilt, X_batchFilt] = ode45(@(t,X)orbitEOM_MuJ2(t,X,pConst.mu, pConst.J2, pConst.Ri), tspan, X0_batch, opt);

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

    %% Create state error plots
boundLevel = 3; % Plot +/- boundLevel*sigma around state errors
titleText = "Batch Filter Estimated State Error (X_{filt} - X_{ref})";
xLabel = "Time [sec]";
yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
          "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]"];

fig_BatchError = plotStateError(t_batchFilt, stateError_batch, t_sigma, sigma_batch, boundLevel, titleText, xLabel, yLabel);

    %% Assign output
batchRun = struct("batchOut", batchOut, "t_batchFilt", t_batchFilt, "X_batchFilt", X_batchFilt, ...
                  "RMS_postfit_batch", RMS_postfit_batch, "RMS_state_comp_batch", RMS_state_comp_batch, ...
                  "RMS_state_full_batch", RMS_state_full_batch, "fig_BatchPreRes", fig_BatchPreRes, ...
                  "fig_BatchPostRes", fig_BatchPostRes, "fig_BatchError", fig_BatchError);

end