function CFARun = runCFA(X0, x0, Pxx0, Pcc0, S0, pConst, stations, X_ref, t_ref, plot)
% Function that runs a CFA on the given data for a stat OD problem. 
%   - Inputs: 
%       - X0: Initial cartesian state from orbital elements, organized as 
%             [X0; Y0; Z0; Xdot0; Ydot0; Zdot0]
%       - x0: Initial state deviation estimate for the filter to use,
%             organized as
%             [x0; y0; z0; xDot0; yDot0; zDot0]
%       - Pxx0: Initial state covariance matrix for the filter to use
%       - Pcc0: Initial consider parameter covariance matrix for the filter
%               to use
%       - S0: Initial sensitivity matrix for the filter to use
%       - pConst: Planetary constants structure as defined by
%                 getPlanetConst.m
%       - stations: Stations structure as defined by makeSations.m
%       - X_ref: Reference orbit from truth data
%       - t_ref: Reference time vector from truth data
%       - numIter: Max number of times to iterate the filter
%       - plot: Boolean indicating whether to plot results or not
%   - Outputs:
%       - CFARun: CFA results structure organized as follows:
%           - CFAOut: Filter output structure as defined in CFA.m
%           - t_CFA: Extracted CFA time vector
%           - X_CFA: Extracted CFA state estimates, organized as follows:
%                    [ 
%                         [X; Y; Z; Xdot; Ydot; Zdot], 
%                         [X; Y; Z; Xdot; Ydot; Zdot], ...
%                    ]
%           - Xc_CFA: Extracted CFA state estimates, including consider 
%                     parameter effects, organized as follows:
%                    [ 
%                         [Xc; Yc; Zc; Xcdot; Ycdot; Zcdot], 
%                         [Xc; Yc; Zc; Xcdot; Ycdot; Zcdot], ...
%                    ]
%           - RMS_prefit_CFA: Prefit residual RMS errors for each run 
%                             of the CFA
%           - RMS_postfit_CFA: Postfit residual RMS errors for each run 
%                              of the CFA
%           - RMS_postfit_c_CFA: Postfit residual RMS errors for each run 
%                                of the CFA, including consider parameter 
%                                effects
%           - RMS_state_comp_CFA: Component-wise state RMS errors, one row 
%                                 per component
%           - RMS_state_full_CFA: Full state RMS error
%           - fig_CFAPreRes: CFA prefit residual plot figure handle
%           - fig_CFAPostRes: CFA postfit residual plot figure handle
%           - fig_CFAPTrace: CFA final P trace plot figure handle
%           - fig_CovEllipsoids: Covariance Ellipsoids plot figure handles
%           - fig_CFAError: CFA state error plot figure handle
%           - fig_CFAError_c: CFA state error plot figure handle, including
%                             consider parameter effects
%
%   By: Ian Faber, 04/07/2025
%
    
    %% Initialize CFA
X0_CFA = X0;
x0_CFA = x0;
Pxx0_CFA = Pxx0;
Pcc0_CFA = Pcc0;
S0_CFA = S0;

fprintf("\n\tRunning CFA:\n")

if any(~plot)
    fig_CFAPreRes = [];
    fig_CFAPostRes = [];
    fig_CFAPTrace = [];
    fig_CovEllipsoids = [];
    fig_CFAError = [];
end

    %% Iterate CFA until residual RMS doesn't change
RMS_prefit_CFA = 1e99;
RMS_postfit_CFA = 1e99; % Start RMS with a bogus value
CFATolerance = 1e-3; % any ratio less than this is considered converged
maxCFARuns = 1; % Cap number of runs
CFARuns = 0;
% fig_CFAPreRes = [];
% fig_CFAPostRes = [];
k = 2; % Start counter with bogus value
while CFARuns < maxCFARuns
        % Run CFA
    CFAOut = CFA(X0_CFA, stations, pConst, Pxx0_CFA, Pcc0_CFA, x0_CFA, S0_CFA);

        % Extract CFA data
    xEst_CFA = CFAOut.xEst;
    PEst_CFA = CFAOut.PEst;
    PcEst_CFA = CFAOut.PcEst;
    prefit_res_CFA = CFAOut.prefit_res;
    postfit_res_CFA = CFAOut.postfit_res;
    t_CFA = CFAOut.t;
    statVis_CFA = CFAOut.statVis;
    X_CFA = CFAOut.XEst;
    Xc_CFA = CFAOut.XcEst;
    Phi_full_CFA = CFAOut.Phi_total{end};
    
        % Calculate residual RMS errors
    rms = calcResidualRMS(prefit_res_CFA, stations, statVis_CFA, true(1,2));
    RMS_prefit_CFA = [RMS_prefit_CFA; rms];
    
    rms = calcResidualRMS(postfit_res_CFA, stations, statVis_CFA, true(1,2));
    RMS_postfit_CFA = [RMS_postfit_CFA; rms];

        % Plot residuals
    % titleText = sprintf("CFA Pre-Fit Residuals - Run %.0f", k-1); 
    % xLabel = "Time [sec]"; 
    % yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
    % colors = ['b', 'r'];
    % fig_CFAPreRes = [fig_CFAPreRes; plotResiduals(t_CFA, prefit_res_CFA, titleText, xLabel, yLabel, colors)];
    % 
    % titleText = sprintf("CFA Post-Fit Residuals - Run %.0f", k-1); 
    % xLabel = "Time [sec]"; 
    % yLabel = ["Range Residuals [m]", "Range-Rate Residuals [m/s]"];
    % colors = ['b', 'r'];
    % fig_CFAPostRes = [fig_CFAPostRes; plotResiduals(t_CFA, postfit_res_CFA, titleText, xLabel, yLabel, colors)];

        % Determine if another run is needed via percent change
    if (abs((RMS_postfit_CFA(k) - RMS_postfit_CFA(k-1))/RMS_postfit_CFA(k-1)) > CFATolerance)
            % Update initial state and perturbations for next run
        x0_CFA = (Phi_full_CFA^-1)*xEst_CFA(:,end);
        X0_CFA = X0_CFA + x0_CFA;

            % Update counters
        k = k + 1;
        CFARuns = CFARuns + 1;

            % Write to console
        if CFARuns < maxCFARuns
            fprintf("Prefit RMS: %.4f, Postfit RMS: %.4f. Iterating CFA. Runs so far: %.0f\n", RMS_prefit_CFA(k-1), RMS_postfit_CFA(k-1), CFARuns)
        else
            fprintf("Prefit RMS: %.4f, Postfit RMS: %.4f. Hit max CFA iterations. Runs so far: %.0f\n", RMS_prefit_CFA(k-1), RMS_postfit_CFA(k-1), CFARuns)
        end
    else
        break;
    end
end
RMS_prefit_CFA = RMS_prefit_CFA(2:end); % Get rid of bogus starting RMS value
RMS_postfit_CFA = RMS_postfit_CFA(2:end); % Get rid of bogus starting RMS value

% if CFARuns < maxCFARuns
fprintf("Final prefit RMS: %.4f\n", RMS_prefit_CFA(end))
fprintf("Final postfit RMS: %.4f\n", RMS_postfit_CFA(end))
% else
%     fprintf("Final prefit RMS: %.4f\n", RMS_prefit_CFA(end))
%     fprintf("Final postfit RMS: %.4f\n", RMS_postfit_CFA(end))
% end

    %% Plot residuals and covariance trace
if plot(1)
            % Residuals
    titleText = sprintf("CFA Pre-Fit Residuals - Run %.0f", CFARuns); 
    xLabel = "Time [sec]"; 
    yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
    colors = ['b', 'r'];
    fig_CFAPreRes = plotResiduals(t_CFA, prefit_res_CFA, titleText, xLabel, yLabel, colors);
end

if plot(2)
    titleText = sprintf("CFA Post-Fit Residuals - Run %.0f", CFARuns); 
    xLabel = "Time [sec]"; 
    yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
    colors = ['b', 'r'];
    fig_CFAPostRes = plotResiduals(t_CFA, postfit_res_CFA, titleText, xLabel, yLabel, colors);
end

if plot(3)
            % Covariance trace
    titleText = sprintf("CFA R Covariance Trace - Run %.0f", CFARuns); 
    xLabel = "Time [sec]"; 
    yLabel = "trace(P)";
    colors = ['b', 'r'];
    elements = 1:3; % Only plot trace of position
    fig_CFAPTrace = plotPTrace(t_CFA, CFAOut.PEst, elements, titleText, xLabel, yLabel, colors);
end

if plot(4)
    titleText = sprintf("CFA V Covariance Trace - Run %.0f", CFARuns); 
    xLabel = "Time [sec]"; 
    yLabel = "trace(P)";
    colors = ['b', 'r'];
    elements = 4:6; % Only plot trace of velocity
    fig_CFAPTrace = [fig_CFAPTrace; plotPTrace(t_CFA, CFAOut.PEst, elements, titleText, xLabel, yLabel, colors)];
end
        %% Plot covariance ellipsoids
    fig_CovEllipsoids = [];
    P_end = PEst_CFA{end};
if plot(5)
    elements = 1:3;
    P_pos = P_end(elements, elements); mu = X_CFA(elements,end);
    titleText = sprintf("Final CFA Position Covariance Ellipsoid, t = %.3f sec\n\\mu = [%.3e, %.3e, %.3e]^T km\n\\sigma_X = %.3e km, \\sigma_Y = %.3e km, \\sigma_Z = %.3e km", t_CFA(end), mu, sqrt(P_pos(1,1)), sqrt(P_pos(2,2)), sqrt(P_pos(3,3)));
    xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; cBarText = "||R - \mu||_2 [km]";
    fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_pos, mu, titleText, xLabel, yLabel, zLabel, cBarText)];
end

if plot(6)
    elements = 4:6;
    P_vel = P_end(elements, elements); mu = X_CFA(elements,end);
    titleText = sprintf("Final CFA Velocity Covariance Ellipsoid, t = %.3f sec\n\\mu = [%.3e, %.3e, %.3e]^T km/s\n\\sigma_{Xdot} = %.3e km/s, \\sigma_{Ydot} = %.3e km/s, \\sigma_{Zdot} = %.3e km/s", t_CFA(end), mu, sqrt(P_vel(1,1)), sqrt(P_vel(2,2)), sqrt(P_vel(3,3)));
    xLabel = "Xdot [km/s]"; yLabel = "Ydot [km/s]"; zLabel = "Zdot [km/s]"; cBarText = "||V - \mu||_2 [km/s]";
    fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_vel, mu, titleText, xLabel, yLabel, zLabel, cBarText)];
end

    %% Extract nominal state from CFA deviation estimates
X_ref_CFA = [];
for k = 1:length(t_CFA)
    X_ref_CFA = [X_ref_CFA; X_ref(t_ref == t_CFA(k),:)];
end

    %% Calculate state error and uncertainty
stateError_CFA = X_CFA' - X_ref_CFA;
stateError_CFA_c = Xc_CFA' - X_ref_CFA;

sigma_CFA = [];
sigma_c_CFA = [];
for k = 1:length(PEst_CFA)
    P = PEst_CFA{k};
    Pc = PcEst_CFA{k};
        
    sigPart = [];
    sigPart_c = [];
    for kk = 1:size(P,1)
        sigPart = [sigPart, sqrt(P(kk,kk))];
        sigPart_c = [sigPart_c, sqrt(Pc(kk,kk))];
    end

    sigma_CFA = [sigma_CFA; sigPart];
    sigma_c_CFA = [sigma_c_CFA; sigPart_c];
end

    %% Find state RMS error: component-wise and state-wise
[RMS_state_comp_CFA, RMS_state_full_CFA] = calcStateErrorRMS(stateError_CFA);

fprintf("State Error Component-wise RMS, no consider effects: [%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]\n", RMS_state_comp_CFA);
fprintf("State Error 3D RMS, no consider effects: %.3e\n", RMS_state_full_CFA);

[RMS_state_comp_CFA_c, RMS_state_full_CFA_c] = calcStateErrorRMS(stateError_CFA_c);

fprintf("State Error Component-wise RMS, including consider effects: [%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]\n", RMS_state_comp_CFA_c);
fprintf("State Error 3D RMS, including consider effects: %.3e\n", RMS_state_full_CFA_c);

    %% Create state error plots
idx = 395;
if plot(7)
    boundLevel = 2; % Plot +/- boundLevel*sigma around state errors
    titleText = "CFA Estimated State Error, ignorning consider effects (X_{filt} - X_{ref})";
    xLabel = "Time [sec]";
    yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
              "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]"];
    
    % fig_CFAError = plotStateError(t_CFA, stateError_CFA, [], [], boundLevel, titleText, xLabel, yLabel);
    fig_CFAError = plotStateError(t_CFA(idx:end), stateError_CFA(idx:end,:), t_CFA(idx:end), sigma_CFA(idx:end,:), boundLevel, titleText, xLabel, yLabel);
end

if plot(8)
    boundLevel = 2; % Plot +/- boundLevel*sigma around state errors
    titleText = "CFA Estimated State Error, including consider effects (X_{filt,c} - X_{ref})";
    xLabel = "Time [sec]";
    yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
              "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]"];
    
    % fig_CFAError = plotStateError(t_CFA, stateError_CFA, [], [], boundLevel, titleText, xLabel, yLabel);
    fig_CFAError_c = plotStateError(t_CFA(idx:end), stateError_CFA(idx:end,:), t_CFA(idx:end), sigma_c_CFA(idx:end,:), boundLevel, titleText, xLabel, yLabel);
end

    %% Assign output
CFARun = struct("CFAOut", CFAOut, "t_CFA", t_CFA, "X_CFA", X_CFA, "Xc_CFA", Xc_CFA, ...
                "RMS_prefit_CFA", RMS_prefit_CFA, "RMS_postfit_CFA", RMS_postfit_CFA, ...
                "RMS_state_comp_CFA", RMS_state_comp_CFA, "RMS_state_full_CFA", RMS_state_full_CFA, ...
                "fig_CFAPreRes", fig_CFAPreRes, "fig_CFAPostRes", fig_CFAPostRes, ...
                "fig_CFAPTrace", fig_CFAPTrace, "fig_CovEllipsoids", fig_CovEllipsoids, ...
                "fig_CFAError", fig_CFAError, "fig_CFAError_c", fig_CFAError_c);

end