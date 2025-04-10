function UKFRun = runUKF(X0, P0, Q0, pConst, stations, X_ref, t_ref, alpha, beta, includeJ3, plot)
% Function that runs a UKF on the given data for a stat OD problem. 
%   - Inputs: 
%       - X0: Initial cartesian state from orbital elements, organized as 
%             [X0; Y0; Z0; Xdot0; Ydot0; Zdot0]
%       - P0: Initial state covariance matrix for the filter to use
%       - Q0: Initial process noise covariance matrix
%       - pConst: Planetary constants structure as defined by
%                 getPlanetConst.m
%       - stations: Stations structure as defined by makeSations.m
%       - X_ref: Reference orbit from truth data
%       - t_ref: Reference time vector from truth data
%       - alpha: UKF sigma point spacing variable from [1e-4, 1]
%       - beta: UKF probability distribution variable, generally 2 for 
%               Gaussian probability distributions
%       - includeJ3: Boolean indicating whether the filter dynamics should
%                    include J3 in addition to mu and J2
%       - plot: Boolean indicating whether to plot results or not
%   - Outputs:
%       - UKFRun: UKF results structure organized as follows:
%           - UKFOut: Filter output structure as defined in UKF.m
%           - t_UKF: Extracted UKF time vector
%           - X_UKF: Extracted UKF state estimates, organized as follows:
%                    [ 
%                         [X; Y; Z; Xdot; Ydot; Zdot], 
%                         [X; Y; Z; Xdot; Ydot; Zdot], ...
%                    ]
%           - RMS_prefit_UKF: Prefit residual RMS error
%           - RMS_postfit_UKF: Postfit residual RMS error
%           - RMS_state_comp_UKF: Component-wise state RMS errors, one row 
%                                 per component
%           - RMS_state_full_UKF: Full state RMS error
%           - fig_UKFPreRes: UKF prefit residual plot figure handle
%           - fig_UKFPostRes: UKF postfit residual plot figure handle
%           - fig_UKFPTrace: UKF final P trace plot figure handle
%           - fig_CovEllipsoids: Covariance Ellipsoids plot figure handles
%           - fig_UKFError: UKF state error plot figure handle
%
%   By: Ian Faber, 02/03/2025
%

    %% Initialize UKF
X0_UKF = X0;
P0_UKF = P0;
Q0_UKF = Q0;

fprintf("\n\tRunning UKF:\n")

if any(~plot)
    fig_UKFPreRes = [];
    fig_UKFPostRes = [];
    fig_UKFPTrace = [];
    fig_CovEllipsoids = [];
    fig_UKFError = [];
end

    %% Run UKF
UKFOut = UKF(stations, pConst, X0_UKF, P0_UKF, Q0_UKF, alpha, beta, includeJ3);

    %% Extract UKF data
XEst_UKF = UKFOut.XEst;
PEst_UKF = UKFOut.PEst;
prefit_res_UKF = UKFOut.prefit_res;
postfit_res_UKF = UKFOut.postfit_res;
t_UKF = UKFOut.t;
statVis_UKF = UKFOut.statVis;

    %% Calculate residual RMS errors
RMS_prefit_UKF = calcResidualRMS(prefit_res_UKF, stations, statVis_UKF, true(1,2));
RMS_postfit_UKF = calcResidualRMS(postfit_res_UKF, stations, statVis_UKF, true(1,2));

fprintf("Prefit RMS: %.4f\n", RMS_postfit_UKF);
fprintf("Postfit RMS: %.4f\n", RMS_postfit_UKF);

    %% Plot residuals and Covariance Trace
if plot(1)
            % Residuals
    titleText = sprintf("UKF Pre-Fit Residuals"); 
    xLabel = "Time [sec]"; 
    yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
    colors = ['b', 'r'];
    fig_UKFPreRes = plotResiduals(t_UKF, real(prefit_res_UKF), titleText, xLabel, yLabel, colors);
end

if plot(2)
    titleText = sprintf("UKF Post-Fit Residuals"); 
    xLabel = "Time [sec]"; 
    yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
    colors = ['b', 'r'];
    fig_UKFPostRes = plotResiduals(t_UKF, real(postfit_res_UKF), titleText, xLabel, yLabel, colors);
end

if plot(3)
            % Covariance trace
    titleText = sprintf("UKF R Covariance Trace"); 
    xLabel = "Time [sec]"; 
    yLabel = "trace(P)";
    colors = ['b', 'r'];
    elements = 1:3; % Only plot trace of position
    fig_UKFPTrace = plotPTrace(t_UKF, UKFOut.PEst, elements, titleText, xLabel, yLabel, colors);
end

if plot(4)
    titleText = sprintf("UKF V Covariance Trace"); 
    xLabel = "Time [sec]"; 
    yLabel = "trace(P)";
    colors = ['b', 'r'];
    elements = 4:6; % Only plot trace of velocity
    fig_UKFPTrace = [fig_UKFPTrace; plotPTrace(t_UKF, UKFOut.PEst, elements, titleText, xLabel, yLabel, colors)];
end

        %% Plot covariance ellipsoids
    fig_CovEllipsoids = [];
    P_end = PEst_UKF{end};

if plot(5)
    elements = 1:3;
    P_pos = P_end(elements, elements); mu = XEst_UKF(elements,end);
    titleText = sprintf("Final UKF Position Covariance Ellipsoid\n\\mu = [%.3e, %.3e, %.3e]^T km\n\\sigma_X = %.3e km, \\sigma_Y = %.3e km, \\sigma_Z = %.3e km", mu, sqrt(P_pos(1,1)), sqrt(P_pos(2,2)), sqrt(P_pos(3,3)));
    xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; cBarText = "||R - \mu||_2 [km]";
    fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_pos, mu, titleText, xLabel, yLabel, zLabel, cBarText)];
end

if plot(6)
    elements = 4:6;
    P_vel = P_end(elements, elements); mu = XEst_UKF(elements,end);
    titleText = sprintf("Final UKF Velocity Covariance Ellipsoid\n\\mu = [%.3e, %.3e, %.3e]^T km/s\n\\sigma_{Xdot} = %.3e km/s, \\sigma_{Ydot} = %.3e km/s, \\sigma_{Zdot} = %.3e km/s", mu, sqrt(P_vel(1,1)), sqrt(P_vel(2,2)), sqrt(P_vel(3,3)));
    xLabel = "Xdot [km/s]"; yLabel = "Ydot [km/s]"; zLabel = "Zdot [km/s]"; cBarText = "||V - \mu||_2 [km/s]";
    fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_vel, mu, titleText, xLabel, yLabel, zLabel, cBarText)];
end

    %% Extract nominal state from UKF deviation estimates
X_ref_UKF = [];
for k = 1:length(t_UKF)
    X_ref_UKF = [X_ref_UKF; X_ref(t_ref == t_UKF(k),:)];
end

    %% Calculate state error and uncertainty
stateError_UKF = XEst_UKF' - X_ref_UKF;

sigma_UKF = [];
for k = 1:length(PEst_UKF)
    P = PEst_UKF{k};

    sigPart = [];
    for kk = 1:size(P,1)
        sigPart = [sigPart, sqrt(P(kk,kk))];
    end

    sigma_UKF = [sigma_UKF; sigPart];
end

    %% Find state RMS error: component-wise and state-wise
[RMS_state_comp_UKF, RMS_state_full_UKF] = calcStateErrorRMS(stateError_UKF);

    %% Create state error plots
if plot(7) || plot(8)
    boundLevel = 3; % Plot +/- boundLevel*sigma around state errors
    titleText = "UKF Estimated State Error (X_{filt} - X_{ref})";
    xLabel = "Time [sec]";
    yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
              "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]"];
    if plot(7)
        fig_UKFError = plotStateError(t_UKF, stateError_UKF, [], [], boundLevel, titleText, xLabel, yLabel);
    end
    
    if plot(8)
        fig_UKFError = plotStateError(t_UKF, stateError_UKF, t_UKF, sigma_UKF, boundLevel, titleText, xLabel, yLabel);
    end
end
    %% Assign output
UKFRun = struct("UKFOut", UKFOut, "t_UKF", t_UKF, "X_UKF", XEst_UKF, ...
                "RMS_prefit_UKF", RMS_prefit_UKF, "RMS_postfit_UKF", RMS_postfit_UKF, ...
                "RMS_state_comp_UKF", RMS_state_comp_UKF, "RMS_state_full_UKF", RMS_state_full_UKF, ...
                "fig_UKFPreRes", fig_UKFPreRes, "fig_UKFPostRes", fig_UKFPostRes, ...
                "fig_UKFPTrace", fig_UKFPTrace, "fig_CovEllipsoids", fig_CovEllipsoids, ...
                "fig_UKFError", fig_UKFError);

end