function SRIFRun = runSRIF(X0, x0, P0, Q0, uBar, forceUpperTriangular, pConst, stations, X_ref, t_ref, numIter, plot)
% Function that runs an SRIF on the given data for a stat OD problem. 
% Iterates SRIF runs until convergence, or the maximum number of iterations 
% is met.
%   - Inputs: 
%       - X0: Initial cartesian state from orbital elements, organized as 
%             [X0; Y0; Z0; Xdot0; Ydot0; Zdot0]
%       - x0: Initial state deviation estimate for the filter to use,
%             organized as
%             [x0; y0; z0; xDot0; yDot0; zDot0]
%       - P0: Initial state covariance matrix for the filter to use
%       - Q0: Initial process noise covariance matrix
%       - uBar: Mean noise vector, generally zeros(3,1) but not always!
%       - forceUpperTriangular: Boolean indicating whether the time updated
%                               R matrix is forced to be upper triangular
%                               or not. Nominally, this should be true.
%       - pConst: Planetary constants structure as defined by
%                 getPlanetConst.m
%       - stations: Stations structure as defined by makeSations.m
%       - X_ref: Reference orbit from truth data
%       - t_ref: Reference time vector from truth data
%       - numIter: Max number of times to iterate the filter
%       - plot: Boolean indicating whether to plot results or not
%   - Outputs:
%       - SRIFRun: SRIF results structure organized as follows:
%           - SRIFOut: Filter output structure as defined in SRIF.m
%           - t_SRIF: Extracted SRIF time vector
%           - X_SRIF: Extracted SRIF state estimates, organized as follows:
%                     [ 
%                          [X; Y; Z; Xdot; Ydot; Zdot], 
%                          [X; Y; Z; Xdot; Ydot; Zdot], ...
%                     ]
%           - RMS_prefit_SRIF_whitened: Whitened prefit residual RMS errors 
%                                       for each run of the SRIF, one run's
%                                       error per row
%           - RMS_postfit_SRIF_whitened: Whitened postfit residual RMS 
%                                        errors for each run of the SRIF, 
%                                        one run's error per row
%           - RMS_prefit_SRIF: Prefit residual RMS errors for each run 
%                              of the SRIF, one run's error per row
%           - RMS_postfit_SRIF: Postfit residual RMS errors for each run 
%                               of the SRIF, one run's error per row
%           - RMS_state_comp_SRIF: Component-wise state RMS errors after 
%                                  the filter has converged, one row per 
%                                  component
%           - RMS_state_full_SRIF: Full state RMS error after the filter
%                                  has converged
%           - fig_SRIFPreRes: SRIF prefit residual plot figure handle
%           - fig_SRIFPostRes: SRIF postfit residual plot figure handle
%           - fig_SRIFPTrace: SRIF final P trace plot figure handle
%           - fig_CovEllipsoids: Covariance Ellipsoids plot figure handles
%           - fig_SRIFError: SRIF state error plot figure handle
%
%   By: Ian Faber, 02/03/2025
%
    
    %% Initialize SRIF
X0_SRIF = X0;
x0_SRIF = x0;
P0_SRIF = P0;
Q0_SRIF = Q0;

fprintf("\n\tRunning SRIF:\n")

if ~plot
    fig_SRIFPreRes = [];
    fig_SRIFPostRes = [];
    fig_SRIFPTrace = [];
    fig_CovEllipsoids = [];
    fig_SRIFError = [];
end

    %% Iterate SRIF until residual RMS doesn't change
RMS_prefit_SRIF_whitened = 1e99;
RMS_postfit_SRIF_whitened = 1e99;
RMS_prefit_SRIF = 1e99;
RMS_postfit_SRIF = 1e99; % Start RMS with a bogus value
SRIFTolerance = 1e-3; % any ratio less than this is considered converged
maxSRIFRuns = numIter; % Cap number of runs
SRIFRuns = 0;
% fig_SRIFPreRes = [];
% fig_SRIFPostRes = [];
k = 2; % Start counter with bogus value
while SRIFRuns < maxSRIFRuns
        % Run SRIF
    SRIFOut = SRIF(X0_SRIF, stations, pConst, P0_SRIF, x0_SRIF, Q0_SRIF, uBar, forceUpperTriangular);

        % Extract SRIF data
    xEst_SRIF = SRIFOut.xEst;
    PEst_SRIF = SRIFOut.PEst;
    prefit_res_SRIF_whitened = SRIFOut.prefit_res_whitened;
    postfit_res_SRIF_whitened = SRIFOut.postfit_res_whitened;
    prefit_res_SRIF = SRIFOut.prefit_res;
    postfit_res_SRIF = SRIFOut.postfit_res;
    t_SRIF = SRIFOut.t;
    statVis_SRIF = SRIFOut.statVis;
    X_SRIF = SRIFOut.XEst;
    Phi_full_SRIF = SRIFOut.Phi_total{end};

    stations_whitened = stations;
    for kk = 1:length(stations_whitened)
        for idx = 1:length(stations_whitened(k).R)
            stations_whitened(kk).R{idx} = eye(size(stations_whitened(kk).R{idx})); %chol(stations_whitened(k).R{kk},'lower'); % Whitened covariance
        end
    end
    
        % Calculate residual RMS errors
    rms = calcResidualRMS(prefit_res_SRIF_whitened, stations_whitened, statVis_SRIF, true(1,2));
    RMS_prefit_SRIF_whitened = [RMS_prefit_SRIF_whitened; rms];
    
    rms = calcResidualRMS(postfit_res_SRIF_whitened, stations_whitened, statVis_SRIF, true(1,2));
    RMS_postfit_SRIF_whitened = [RMS_postfit_SRIF_whitened; rms];

    rms = calcResidualRMS(prefit_res_SRIF, stations, statVis_SRIF, true(1,2));
    RMS_prefit_SRIF = [RMS_prefit_SRIF; rms];
    
    rms = calcResidualRMS(postfit_res_SRIF, stations, statVis_SRIF, true(1,2));
    RMS_postfit_SRIF = [RMS_postfit_SRIF; rms];

        % Plot residuals
    % titleText = sprintf("SRIF Pre-Fit Residuals - Run %.0f", k-1); 
    % xLabel = "Time [sec]"; 
    % yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
    % colors = ['b', 'r'];
    % fig_SRIFPreRes = [fig_SRIFPreRes; plotResiduals(t_SRIF, prefit_res_SRIF, titleText, xLabel, yLabel, colors)];
    % 
    % titleText = sprintf("SRIF Post-Fit Residuals - Run %.0f", k-1); 
    % xLabel = "Time [sec]"; 
    % yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
    % colors = ['b', 'r'];
    % fig_SRIFPostRes = [fig_SRIFPostRes; plotResiduals(t_SRIF, postfit_res_SRIF, titleText, xLabel, yLabel, colors)];

        % Determine if another run is needed via percent change
    if (abs((RMS_postfit_SRIF(k) - RMS_postfit_SRIF(k-1))/RMS_postfit_SRIF(k-1)) > SRIFTolerance)
            % Update initial state and perturbations for next run
        x0_SRIF = (Phi_full_SRIF^-1)*xEst_SRIF(:,end);
        X0_SRIF = X0_SRIF + x0_SRIF;

            % Update counters
        k = k + 1;
        SRIFRuns = SRIFRuns + 1;

            % Write to console
        if SRIFRuns < maxSRIFRuns
            fprintf("Prefit RMS: %.4f, Postfit RMS: %.4f. Iterating SRIF. Runs so far: %.0f\n", RMS_prefit_SRIF(k-1), RMS_postfit_SRIF(k-1), SRIFRuns)
        else
            fprintf("Prefit RMS: %.4f, Postfit RMS: %.4f. Hit max SRIF iterations. Runs so far: %.0f\n", RMS_prefit_SRIF(k-1), RMS_postfit_SRIF(k-1), SRIFRuns)
        end
    else
        break;
    end
end
RMS_prefit_SRIF_whitened = RMS_prefit_SRIF_whitened(2:end); % Get rid of bogus starting RMS value
RMS_postfit_SRIF_whitened = RMS_postfit_SRIF_whitened(2:end); % Get rid of bogus starting RMS value
RMS_prefit_SRIF = RMS_prefit_SRIF(2:end); % Get rid of bogus starting RMS value
RMS_postfit_SRIF = RMS_postfit_SRIF(2:end); % Get rid of bogus starting RMS value

if SRIFRuns < maxSRIFRuns
    fprintf("Final whitened prefit RMS: %.4f. Converged after %.0f runs\n", RMS_prefit_SRIF_whitened(end), SRIFRuns)
    fprintf("Final whitened postfit RMS: %.4f. Converged after %.0f runs\n", RMS_postfit_SRIF_whitened(end), SRIFRuns)
    fprintf("Final prefit RMS: %.4f. Converged after %.0f runs\n", RMS_prefit_SRIF(end), SRIFRuns)
    fprintf("Final postfit RMS: %.4f. Converged after %.0f runs\n", RMS_postfit_SRIF(end), SRIFRuns)
else
    fprintf("Final whitened prefit RMS: %.4f. Hit maximum number of %.0f runs\n", RMS_prefit_SRIF_whitened(end), maxSRIFRuns)
    fprintf("Final whitened postfit RMS: %.4f. Hit maximum number of %.0f runs\n", RMS_postfit_SRIF_whitened(end), maxSRIFRuns)
    fprintf("Final prefit RMS: %.4f. Hit maximum number of %.0f runs\n", RMS_prefit_SRIF(end), maxSRIFRuns)
    fprintf("Final postfit RMS: %.4f. Hit maximum number of %.0f runs\n", RMS_postfit_SRIF(end), maxSRIFRuns)
end

    %% Plot residuals and covariance trace
if plot
            % Residuals
    titleText = sprintf("SRIF Whitened Pre-Fit Residuals - Run %.0f", SRIFRuns); 
    xLabel = "Time [sec]"; 
    yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
    colors = ['b', 'r'];
    fig_SRIFPreRes_Whitened = plotResiduals(t_SRIF, prefit_res_SRIF_whitened, titleText, xLabel, yLabel, colors);
    
    titleText = sprintf("SRIF Whitened Post-Fit Residuals - Run %.0f", SRIFRuns); 
    xLabel = "Time [sec]"; 
    yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
    colors = ['b', 'r'];
    fig_SRIFPostRes_Whitened = plotResiduals(t_SRIF, postfit_res_SRIF_whitened, titleText, xLabel, yLabel, colors);

    titleText = sprintf("SRIF Pre-Fit Residuals - Run %.0f", SRIFRuns); 
    xLabel = "Time [sec]"; 
    yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
    colors = ['b', 'r'];
    fig_SRIFPreRes = plotResiduals(t_SRIF, prefit_res_SRIF, titleText, xLabel, yLabel, colors);
    
    titleText = sprintf("SRIF Post-Fit Residuals - Run %.0f", SRIFRuns); 
    xLabel = "Time [sec]"; 
    yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
    colors = ['b', 'r'];
    fig_SRIFPostRes = plotResiduals(t_SRIF, postfit_res_SRIF, titleText, xLabel, yLabel, colors);
            
            % Covariance trace
    titleText = sprintf("SRIF R Covariance Trace - Run %.0f", SRIFRuns); 
    xLabel = "Time [sec]"; 
    yLabel = "trace(P)";
    colors = ['b', 'r'];
    elements = 1:3; % Only plot trace of position
    fig_SRIFPTrace = plotPTrace(t_SRIF, SRIFOut.PEst, elements, titleText, xLabel, yLabel, colors);
    
    titleText = sprintf("SRIF V Covariance Trace - Run %.0f", SRIFRuns); 
    xLabel = "Time [sec]"; 
    yLabel = "trace(P)";
    colors = ['b', 'r'];
    elements = 4:6; % Only plot trace of velocity
    fig_SRIFPTrace = [fig_SRIFPTrace; plotPTrace(t_SRIF, SRIFOut.PEst, elements, titleText, xLabel, yLabel, colors)];
    
        %% Plot covariance ellipsoids
    fig_CovEllipsoids = [];
    P_end = PEst_SRIF{end};
    
    elements = 1:3;
    P_pos = P_end(elements, elements); mu = X_SRIF(elements,end);
    titleText = sprintf("Final SRIF Position Covariance Ellipsoid, t = %.3f sec\n\\mu = [%.3e, %.3e, %.3e]^T km\n\\sigma_X = %.3e km, \\sigma_Y = %.3e km, \\sigma_Z = %.3e km", t_SRIF(end), mu, sqrt(P_pos(1,1)), sqrt(P_pos(2,2)), sqrt(P_pos(3,3)));
    xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; cBarText = "||R - \mu||_2 [km]";
    fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_pos, mu, titleText, xLabel, yLabel, zLabel, cBarText)];
    
    elements = 4:6;
    P_vel = P_end(elements, elements); mu = X_SRIF(elements,end);
    titleText = sprintf("Final SRIF Velocity Covariance Ellipsoid, t = %.3f sec\n\\mu = [%.3e, %.3e, %.3e]^T km/s\n\\sigma_{Xdot} = %.3e km/s, \\sigma_{Ydot} = %.3e km/s, \\sigma_{Zdot} = %.3e km/s", t_SRIF(end), mu, sqrt(P_vel(1,1)), sqrt(P_vel(2,2)), sqrt(P_vel(3,3)));
    xLabel = "Xdot [km/s]"; yLabel = "Ydot [km/s]"; zLabel = "Zdot [km/s]"; cBarText = "||V - \mu||_2 [km/s]";
    fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_vel, mu, titleText, xLabel, yLabel, zLabel, cBarText)];
end

    %% Extract nominal state from SRIF deviation estimates
X_ref_SRIF = [];
for k = 1:length(t_SRIF)
    X_ref_SRIF = [X_ref_SRIF; X_ref(t_ref == t_SRIF(k),:)];
end

    %% Calculate state error and uncertainty
stateError_SRIF = X_SRIF' - X_ref_SRIF;

sigma_SRIF = [];
for k = 1:length(PEst_SRIF)
    P = PEst_SRIF{k};

    sigPart = [];
    for kk = 1:size(P,1)
        sigPart = [sigPart, sqrt(P(kk,kk))];
    end

    sigma_SRIF = [sigma_SRIF; sigPart];
end

    %% Find state RMS error: component-wise and state-wise
[RMS_state_comp_SRIF, RMS_state_full_SRIF] = calcStateErrorRMS(stateError_SRIF);

    %% Create state error plots
if plot
    boundLevel = 3; % Plot +/- boundLevel*sigma around state errors
    titleText = "SRIF Estimated State Error (X_{filt} - X_{ref})";
    xLabel = "Time [sec]";
    yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
              "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]"];
    
    fig_SRIFError = plotStateError(t_SRIF, stateError_SRIF, [], [], boundLevel, titleText, xLabel, yLabel);
    fig_SRIFError = plotStateError(t_SRIF, stateError_SRIF, t_SRIF, sigma_SRIF, boundLevel, titleText, xLabel, yLabel);
end
    %% Assign output
SRIFRun = struct("SRIFOut", SRIFOut, "t_SRIF", t_SRIF, "X_SRIF", X_SRIF, ...
                "RMS_prefit_SRIF_whitened", RMS_prefit_SRIF_whitened, "RMS_postfit_SRIF_whitened", RMS_postfit_SRIF_whitened, ...
                "RMS_prefit_SRIF", RMS_prefit_SRIF, "RMS_postfit_SRIF", RMS_postfit_SRIF, ...
                "RMS_state_comp_SRIF", RMS_state_comp_SRIF, "RMS_state_full_SRIF", RMS_state_full_SRIF, ...
                "fig_SRIFPreRes_Whitened", fig_SRIFPreRes_Whitened, "fig_SRIFPostRes_Whitened", fig_SRIFPostRes_Whitened, ...
                "fig_SRIFPreRes", fig_SRIFPreRes, "fig_SRIFPostRes", fig_SRIFPostRes, ...
                "fig_SRIFPTrace", fig_SRIFPTrace, "fig_CovEllipsoids", fig_CovEllipsoids, ...
                "fig_SRIFError", fig_SRIFError);

end