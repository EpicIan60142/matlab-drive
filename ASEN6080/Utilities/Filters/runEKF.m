function EKFRun = runEKF(X0, P0, pConst, stations, X_ref, t_ref, t_start)
% Function that runs an LKF on the given data for a stat OD problem. 
%   - Inputs: 
%       - X0: Initial cartesian state from orbital elements, organized as 
%             [X0; Y0; Z0; Xdot0; Ydot0; Zdot0]
%       - P0: Initial state covariance matrix for the filter to use
%       - pConst: Planetary constants structure as defined by
%                 getPlanetConst.m
%       - stations: Stations structure as defined by makeSations.m
%       - X_ref: Reference orbit from truth data
%       - t_ref: Reference time vector from truth data
%       - t_start: Time for the EKF to start processing measurements at,
%                  often tied to the state used to initialize the EKF, i.e.
%                  the time after the last measurement of an LKF.
%   - Outputs:
%       - EKFRun: EKF results structure organized as follows:
%           - EKFOut: Filter output structure as defined in EKF.m
%           - t_EKF: Extracted EKF time vector
%           - X_EKF: Extracted EKF state estimates, organized as follows:
%                    [ 
%                         [X; Y; Z; Xdot; Ydot; Zdot], 
%                         [X; Y; Z; Xdot; Ydot; Zdot], ...
%                    ]
%           - RMS_postfit_EKF: Postfit residual RMS error
%           - RMS_state_comp_EKF: Component-wise state RMS errors, one row 
%                                 per component
%           - RMS_state_full_EKF: Full state RMS error
%           - fig_EKFRes: EKF residual plot figure handle
%           - fig_EKFError: EKF state error plot figure handle
%
%   By: Ian Faber, 02/03/2025
%

    %% Initialize EKF
X0_EKF = X0;
P0_EKF = P0;
t_EKF_start = t_start;

fprintf("\n\tRunning EKF:\n")

    %% Run EKF
EKFOut = EKF(X0_EKF, stations, pConst, P0_EKF, t_EKF_start);

    %% Extract EKF data
xEst_EKF = EKFOut.xEst;
PEst_EKF = EKFOut.PEst;
prefit_res_EKF = EKFOut.prefit_res;
postfit_res_EKF = EKFOut.postfit_res;
t_EKF = EKFOut.t;
statVis_EKF = EKFOut.statVis;
X_EKF = EKFOut.XEst;

    %% Calculate residual RMS errors
RMS_postfit_EKF = calcResidualRMS(postfit_res_EKF, stations, statVis_EKF);

fprintf("Postfit RMS: %.4f\n", RMS_postfit_EKF);

    %% Plot residuals
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

    %% Extract nominal state from LKF deviation estimates
X_ref_EKF = [];
for k = 1:length(t_EKF)
    X_ref_EKF = [X_ref_EKF; X_ref(t_ref == t_EKF(k),:)];
end

    %% Calculate state error and uncertainty
stateError_EKF = X_EKF' - X_ref_EKF;

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

    %% Create state error plots
boundLevel = 3; % Plot +/- boundLevel*sigma around state errors
titleText = "EKF Estimated State Error (X_{filt} - X_{ref})";
xLabel = "Time [sec]";
yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
          "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]"];

fig_EKFError = plotStateError(t_EKF, stateError_EKF, t_EKF, sigma_EKF, boundLevel, titleText, xLabel, yLabel);

    %% Assign output
EKFRun = struct("EKFOut", EKFOut, "t_EKF", t_EKF, "X_EKF", X_EKF, ...
                "RMS_postfit_EKF", RMS_postfit_EKF, "RMS_state_comp_EKF", RMS_state_comp_EKF, ...
                "RMS_state_full_EKF", RMS_state_full_EKF, "fig_EKFRes", fig_EKFPreRes, ...
                "fig_EKFPostRes", fig_EKFPostRes, "fig_EKFError", fig_EKFError);

end