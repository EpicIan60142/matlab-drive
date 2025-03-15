function smoothRun = runSmoother(filterOut, X_ref, t_ref, plot)
% Function that runs a sequential smoother algorithm and plots outputs for
% Stat OD problems
%   - Inputs:
%       - filterOut: Output structure from an LKF as defined by LKF_SNC.m
%       - X_ref: Reference orbit from truth data
%       - t_ref: Reference time vector from truth data
%       - plot: Boolean indicating whether to plot results or not
%   - Outputs:
%       - smoothRun: Smoother run output structure with the following
%                    fields:
%           - smoothOut: Smoother output structure as defined in Smoother.m
%           - t_Smooth: Extracted Smoother time vector
%           - X_Smooth: Extracted Smoother state estimates, organized as 
%                       follows:
%                       [ 
%                             [X; Y; Z; Xdot; Ydot; Zdot], 
%                             [X; Y; Z; Xdot; Ydot; Zdot], ...
%                       ]
%           - RMS_state_comp_Smooth: Component-wise state RMS errors, one 
%                                    row per component
%           - RMS_state_full_Smooth: Full state RMS error
%           - fig_SmoothPTrace: Smoother final P trace plot figure handle
%           - fig_CovEllipsoids: Covariance Ellipsoids plot figure handles
%           - fig_SmoothError: Smoother state error plot figure handle
    %% Initialize outputs if needed
if ~plot
    fig_SmoothPTrace = [];
    fig_CovEllipsoids = [];
    fig_SmoothError = [];
end

fprintf("\n\tRunning Smoother...\n")
    
    %% Run smoother
smoothOut = Smoother(filterOut);

    %% Extract results
t_Smooth = smoothOut.tSmoothed;
X_Smooth = smoothOut.XSmoothed;
P_Smooth = smoothOut.PSmoothed;

    %% Plot covariance results
if plot
        % Covariance trace
    titleText = sprintf("Smoother R Covariance Trace"); 
    xLabel = "Time [sec]"; 
    yLabel = "trace(P)";
    colors = ['b', 'r'];
    elements = 1:3; % Only plot trace of position
    fig_SmoothPTrace = plotPTrace(t_Smooth, smoothOut.PSmoothed, elements, titleText, xLabel, yLabel, colors);
    
    titleText = sprintf("Smoother V Covariance Trace"); 
    xLabel = "Time [sec]"; 
    yLabel = "trace(P)";
    colors = ['b', 'r'];
    elements = 4:6; % Only plot trace of velocity
    fig_SmoothPTrace = [fig_SmoothPTrace; plotPTrace(t_Smooth, smoothOut.PSmoothed, elements, titleText, xLabel, yLabel, colors)];

        % Covariance ellipsoids
    fig_CovEllipsoids = [];
    P_end = P_Smooth{end};
    
    elements = 1:3;
    P_pos = P_end(elements, elements); mu = X_Smooth(elements,end);
    titleText = sprintf("Final Smoother Position Covariance Ellipsoid, t = %.3f sec\n\\mu = [%.3e, %.3e, %.3e]^T km\n\\sigma_X = %.3e km, \\sigma_Y = %.3e km, \\sigma_Z = %.3e km", t_Smooth(end), mu, sqrt(P_pos(1,1)), sqrt(P_pos(2,2)), sqrt(P_pos(3,3)));
    xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; cBarText = "||R - \mu||_2 [km]";
    fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_pos, mu, titleText, xLabel, yLabel, zLabel, cBarText)];
    
    elements = 4:6;
    P_vel = P_end(elements, elements); mu = X_Smooth(elements,end);
    titleText = sprintf("Final Smoother Velocity Covariance Ellipsoid, t = %.3f sec\n\\mu = [%.3e, %.3e, %.3e]^T km/s\n\\sigma_{Xdot} = %.3e km/s, \\sigma_{Ydot} = %.3e km/s, \\sigma_{Zdot} = %.3e km/s", t_Smooth(end), mu, sqrt(P_vel(1,1)), sqrt(P_vel(2,2)), sqrt(P_vel(3,3)));
    xLabel = "Xdot [km/s]"; yLabel = "Ydot [km/s]"; zLabel = "Zdot [km/s]"; cBarText = "||V - \mu||_2 [km/s]";
    fig_CovEllipsoids = [fig_CovEllipsoids; plotCovEllipsoid(P_vel, mu, titleText, xLabel, yLabel, zLabel, cBarText)];
end

    %% Extract nominal state from Smoother deviation estimates
X_ref_Smooth = [];
for k = 1:length(t_Smooth)
    X_ref_Smooth = [X_ref_Smooth; X_ref(t_ref == t_Smooth(k),:)];
end

    %% Calculate state error and uncertainty
stateError_Smooth = X_Smooth' - X_ref_Smooth;

sigma_Smooth = [];
for k = 1:length(P_Smooth)
    P = P_Smooth{k};

    sigPart = [];
    for kk = 1:size(P,1)
        sigPart = [sigPart, sqrt(P(kk,kk))];
    end

    sigma_Smooth = [sigma_Smooth; sigPart];
end

    %% Find state RMS error: component-wise and state-wise
[RMS_state_comp_Smooth, RMS_state_full_Smooth] = calcStateErrorRMS(stateError_Smooth);

    %% Create state error plots
if plot
    boundLevel = 3; % Plot +/- boundLevel*sigma around state errors
    titleText = "Smoother Estimated State Error (X_{smooth} - X_{ref})";
    xLabel = "Time [sec]";
    yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
              "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]"];
    
    fig_SmoothError = plotStateError(t_Smooth, stateError_Smooth, [], [], boundLevel, titleText, xLabel, yLabel);
    fig_SmoothError = plotStateError(t_Smooth, stateError_Smooth, t_Smooth, sigma_Smooth, boundLevel, titleText, xLabel, yLabel);
end
    %% Assign output
smoothRun = struct("smoothOut", smoothOut, "t_Smooth", t_Smooth, ...
                   "X_Smooth", X_Smooth, "RMS_state_comp_Smooth", RMS_state_comp_Smooth, ...
                   "RMS_state_full_Smooth", RMS_state_full_Smooth, ...
                   "fig_SmoothPTrace", fig_SmoothPTrace, "fig_CovEllipsoids", fig_CovEllipsoids, ...
                   "fig_SmoothError", fig_SmoothError);

end

