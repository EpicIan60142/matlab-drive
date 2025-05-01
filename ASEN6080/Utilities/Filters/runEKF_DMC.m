function EKFRun = runEKF_DMC(X0, P0, B, Qu, pConst, stations, X_ref, t_ref, t_start, plot)
% Function that runs an LKF on the given data for a stat OD problem
% including dynamic model compensation
%   - Inputs: 
%       - X0: Initial cartesian state from orbital elements, organized as 
%             [X0; Y0; Z0; Xdot0; Ydot0; Zdot0]
%       - P0: Initial state covariance matrix for the filter to use
%       - B: DMC time constant matrix organized as a 3x3 matrix:
%            B = diag([tau_x^-1, tau_y^-1, tau_z^-1])
%       - Qu: Process noise covariance matrix organized as a 3x3 matrix:
%             Qu = diag([sigma_x^2, sigma_y^2, sigma_z^2])
%       - pConst: Planetary constants structure as defined by
%                 getPlanetConst.m
%       - stations: Stations structure as defined by makeSations.m
%       - X_ref: Reference orbit from truth data
%       - t_ref: Reference time vector from truth data
%       - t_start: Time for the EKF to start processing measurements at,
%                  often tied to the state used to initialize the EKF, i.e.
%                  the time after the last measurement of an LKF.
%       - plot: Boolean indicating whether to plot results or not
%   - Outputs:
%       - EKFRun: EKF results structure organized as follows:
%           - EKFOut: Filter output structure as defined in EKF.m
%           - t_EKF: Extracted EKF time vector
%           - X_EKF: Extracted EKF state estimates, organized as follows:
%                    [ 
%                         [X; Y; Z; Xdot; Ydot; Zdot], 
%                         [X; Y; Z; Xdot; Ydot; Zdot], ...
%                    ]
%           - RMS_prefit_EKF: Prefit residual RMS error
%           - RMS_postfit_EKF: Postfit residual RMS error
%           - RMS_state_comp_EKF: Component-wise state RMS errors, one row 
%                                 per component
%           - RMS_state_full_EKF: Full state RMS error
%           - fig_LKFPreRes: LKF prefit residual plot figure handle
%           - fig_LKFPostRes: LKF postfit residual plot figure handle
%           - fig_LKFPTrace: LKF final P trace plot figure handle
%           - fig_CovEllipsoids: Covariance Ellipsoids plot figure handles
%           - fig_EKFError: EKF state error plot figure handle
%
%   By: Ian Faber, 02/03/2025
%

    %% Initialize EKF
X0_EKF = X0;
P0_EKF = P0;
t_EKF_start = t_start;

fprintf("\n\tRunning EKF:\n")

if ~plot
    fig_EKFPreRes = [];
    fig_EKFPostRes = [];
    fig_EKFPTrace = [];
    fig_CovEllipsoids = [];
    fig_EKFError = [];
end

    %% Run EKF
EKFOut = EKF_DMC(X0_EKF, stations, pConst, P0_EKF, B, Qu, t_EKF_start);

    %% Extract EKF data
xEst_EKF = EKFOut.xEst;
PEst_EKF = EKFOut.PEst;
prefit_res_EKF = EKFOut.prefit_res;
postfit_res_EKF = EKFOut.postfit_res;
t_EKF = EKFOut.t;
statVis_EKF = EKFOut.statVis;
X_EKF = EKFOut.XEst;

    %% Calculate residual RMS errors
RMS_prefit_EKF = calcResidualRMS(prefit_res_EKF, stations, statVis_EKF, true(1,2));
RMS_postfit_EKF = calcResidualRMS(postfit_res_EKF, stations, statVis_EKF, true(1,2));

fprintf("Prefit RMS: %.4f\n", RMS_postfit_EKF);
fprintf("Postfit RMS: %.4f\n", RMS_postfit_EKF);

if plot
    %% Plot residuals and Covariance Trace
        % Residuals
titleText = sprintf("EKF Pre-Fit Residuals"); 
xLabel = "Time [sec]"; 
yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
colors = ['b', 'r'];
fig_EKFPreRes = plotResiduals(t_EKF, prefit_res_EKF, titleText, xLabel, yLabel, colors);

titleText = sprintf("EKF Post-Fit Residuals"); 
    xLabel = "Time [sec]"; 
    yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
    colors = ['b', 'r'];
    fig_EKFPostRes = plotResiduals(t_EKF, postfit_res_EKF, titleText, xLabel, yLabel, colors);
            
            % Covariance trace
    titleText = sprintf("EKF R Covariance Trace"); 
    xLabel = "Time [sec]"; 
    yLabel = "trace(P)";
    colors = ['b', 'r'];
    elements = 1:3; % Only plot trace of position
    fig_EKFPTrace = plotPTrace(t_EKF, EKFOut.PEst, elements, titleText, xLabel, yLabel, colors);
    
    titleText = sprintf("EKF V Covariance Trace"); 
    xLabel = "Time [sec]"; 
    yLabel = "trace(P)";
    colors = ['b', 'r'];
    elements = 4:6; % Only plot trace of velocity
    fig_EKFPTrace = [fig_EKFPTrace; plotPTrace(t_EKF, EKFOut.PEst, elements, titleText, xLabel, yLabel, colors)];
    
        %% Plot covariance ellipsoids
    fig_CovEllipsoids = [];
    P_end = PEst_EKF{end};
    
    elements = 1:3;
    P_pos = P_end(elements, elements); mu = X_EKF(elements,end);
    titleText = sprintf("Final EKF Position Covariance Ellipsoid\n\\mu = [%.3e, %.3e, %.3e]^T km\n\\sigma_X = %.3e km, \\sigma_Y = %.3e km, \\sigma_Z = %.3e m", mu, sqrt(P_pos(1,1)), sqrt(P_pos(2,2)), sqrt(P_pos(3,3)));
    xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; cBarText = "||R - \mu||_2 [km]";
    fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_pos, mu, titleText, xLabel, yLabel, zLabel, cBarText)];
    
    elements = 4:6;
    P_vel = P_end(elements, elements); mu = X_EKF(elements,end);
    titleText = sprintf("Final EKF Velocity Covariance Ellipsoid\n\\mu = [%.3e, %.3e, %.3e]^T km/s\n\\sigma_{Xdot} = %.3e km/s, \\sigma_{Ydot} = %.3e km/s, \\sigma_{Zdot} = %.3e km/s", mu, sqrt(P_vel(1,1)), sqrt(P_vel(2,2)), sqrt(P_vel(3,3)));
    xLabel = "Xdot [km/s]"; yLabel = "Ydot [km/s]"; zLabel = "Zdot [km/s]"; cBarText = "||V - \mu||_2 [km/s]";
    fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_vel, mu, titleText, xLabel, yLabel, zLabel, cBarText)];

end

    %% Extract nominal state from LKF deviation estimates
X_ref_EKF = [];
for k = 1:length(t_EKF)
    X_ref_EKF = [X_ref_EKF; X_ref(t_ref == t_EKF(k),:)];
end

    %% Calculate state error and uncertainty
stateError_EKF = X_EKF(1:6,:)' - X_ref_EKF;

sigma_EKF = [];
for k = 1:length(PEst_EKF)
    P = PEst_EKF{k};

    sigPart = [];
    for kk = 1:size(P,1)
        sigPart = [sigPart, sqrt(P(kk,kk))];
    end

    sigma_EKF = [sigma_EKF; sigPart];
end

    %% Find state RMS error: component-wise and state-wise
[RMS_state_comp_EKF, RMS_state_full_EKF] = calcStateErrorRMS(stateError_EKF);

if plot
        %% Create state error plots
    boundLevel = 3; % Plot +/- boundLevel*sigma around state errors
    titleText = "EKF Estimated State Error (X_{filt} - X_{ref})";
    xLabel = "Time [sec]";
    yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
              "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]"];
    
    fig_EKFError = plotStateError(t_EKF, stateError_EKF, [], [], boundLevel, titleText, xLabel, yLabel);
    fig_EKFError = plotStateError(t_EKF, stateError_EKF, t_EKF, sigma_EKF, boundLevel, titleText, xLabel, yLabel);

end
    %% Assign output
EKFRun = struct("EKFOut", EKFOut, "t_EKF", t_EKF, "X_EKF", X_EKF, ...
                "RMS_prefit_EKF", RMS_prefit_EKF, "RMS_postfit_EKF", RMS_postfit_EKF, ...
                "RMS_state_comp_EKF", RMS_state_comp_EKF, "RMS_state_full_EKF", RMS_state_full_EKF, ...
                "fig_EKFPreRes", fig_EKFPreRes, "fig_EKFPostRes", fig_EKFPostRes, ...
                "fig_EKFPTrace", fig_EKFPTrace, "fig_CovEllipsoids", fig_CovEllipsoids, ...
                "fig_EKFError", fig_EKFError);

end