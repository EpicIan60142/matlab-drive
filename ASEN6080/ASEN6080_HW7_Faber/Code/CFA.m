function filterOut = CFA(Xstar0, stations, pConst, P0, Pcc0, x0, S0, numMeas)
% Function that implements the Consider Filter Algorithm for stat OD 
% problems
%   Inputs:
%       - Xstar0: Initial value of the reference trajectory, organized as
%                 follows:
%                 [X0; Y0; Z0; Xdot0; Ydot0; Zdot0]
%       - stations: Stations struct as defined by makeStations.m. Must have
%                   propagated station states! To propagate states, see
%                   generateTruthData.m.
%       - pConst: Planetary constant structure as formatted by
%                 getPlanetConst.m
%       - P0: Initial state covariance estimate
%       - Pcc0: Initial consider parameter covariance estimate
%       - x0: Initial state deviation estimate
%       - S0: Initial sensitivity analysis
%       - numMeas: Number of measurements to process (optional). If not 
%                  specified, defaults to all measurements
%   Outputs:
%       - filterOut: Filter output structure with the following fields:
%           - xEst: Estimated state deviation at each time processed from 
%                   the station measurements in "stations" (t), organized
%                   as follows:
%                   [xEst_1, xEst_2, ..., x_t], where 
%                   xEst = [x; y; z; xDot; yDot; zDot]
%           - xcEst: Estimated state deviation at each time processed from 
%                    the station measurements in "stations" (t), including 
%                    consider parameter effects, organized as follows:
%                    [xcEst_1, xcEst_2, ..., xc_t], where 
%                    xcEst = [xc; yc; zc; xcDot; ycDot; zcDot]
%           - PEst: Estimated state covariance at each time in t, organized 
%                   as follows:
%                   [{P_1}, {P_2}, ..., {P_t}]
%           - PcEst: Estimated state covariance at each time in t, 
%                    including consider parameter effects, organized as 
%                    follows:
%                    [{Pc_1}, {Pc_2}, ..., {Pc_t}]
%           - PxcEst: Estimated state to consider parameter covariance at 
%                     each time in t, organized as follows:
%                     [{Pxc_1}, {Pxc_2}, ..., {Pxc_t}]
%           - prefit_res: Pre-fit residuals (y_i) at each time in t:
%                         [y_1, y_2, ..., y_t]
%           - postfit_res: Postfit residuals (epsilon = y_i - Htilde_i*x_i)
%                          at each time in t:
%                          [epsilon_1, epsilon_2, ..., epsilon_t]
%           - postfit_res_c: Postfit residuals (epsilon_c = y_i - Htilde_i*x_ci)
%                            at each time in t:
%                            [epsilon_c1, epsilon_c2, ..., epsilon_ct]
%           - t: Measurement time vector for the LKF filter
%           - statVis: Station visibility vector
%           - XEst: Estimated full state at each time in t:
%                   [XEst_1, XEst_2, ..., XEst_t], where
%                   XEst = [X; Y; Z; XDot; YDot; ZDot]
%           - XcEst: Estimated full state at each time in t, including 
%                    consider parameter effects:
%                    [XcEst_1, XcEst_2, ..., XcEst_t], where
%                    XcEst = [Xc; Yc; Zc; XcDot; YcDot; ZcDot]
%           - Phi_total: Integrated STM from t0 to t_i for each time in t:
%                        [{Phi_1}, {Phi_2}, ..., {Phi_t}], where
%                        Phi_t = Phi(t_i, t_0)
%           - Psi: Integrated consider parameter and state mapping matrix
%                  from t_0 to t_i for each time in t:
%                  [{Psi_1}, {Psi_2}, ..., {Psi_t}], where 
%                  Psi_t = Psi(t_i, t_0)
%
%   By: Ian Faber, 04/05/2025
%

%% Initialize settings

    % Format ode45 and sizes
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
n = length(Xstar0);
xEst = [];
xcEst = [];
PEst = [];
PcEst = [];
PxcEst = [];
prefit_res = [];
postfit_res = [];
postfit_res_c = [];
XEst = [];
XcEst = [];
Phi_total = [];
Theta_total = [];
Psi = [];

    %% Process station data into a usable form
[t, Y, R, Xs, vis] = processStations(stations);

    % Specify number of measurements to process
if ~exist("numMeas", 'var')
    numMeas = length(Y);
end

    %% Loop through each observation
t_im1 = t(1);
Xstar_im1 = Xstar0;
x_im1 = x0;
P_im1 = P0;
S_im1 = S0;
Phi_full = eye(n);
Theta_full = zeros(6,1);

[~, XPhi_test] = ode45(@(t,XPhi)STMEOM_J2(t,XPhi,pConst.mu,pConst.J2,pConst.Ri), [t(1), t(end)], [Xstar0; reshape(eye(n),n^2,1)], opt);
Phi_test = reshape(XPhi_test(end,n+1:end), n, n);

[~, XTheta_test] = ode45(@(t,XTheta)STMEOM_CP_J3(t,XTheta,pConst), [t(1), t(end)], [Xstar0; zeros(n,1)], opt);
Theta_test = reshape(XTheta_test(end,n+1:end), n, 1);

for k = 2:numMeas

        % Read next time, measurement, and measurement covariance
    t_i = t(k);
    Y_i = Y{k};
    R_i = R{k};

        % Continue to integrate Phi(t, t0) and Theta(t, t0) for iteration purposes
    XPhi_full = [Xstar_im1; reshape(Phi_full,n^2,1)];
    [~, XPhi_full] = ode45(@(t,XPhi)STMEOM_J2(t,XPhi,pConst.mu, pConst.J2, pConst.Ri), [t_im1 t_i], XPhi_full, opt);
    Phi_full = reshape(XPhi_full(end,n+1:end), n, n);

    XTheta_full = [Xstar_im1; reshape(Theta_full,n,1)];
    [~, XTheta_full] = ode45(@(t,XTheta)STMEOM_CP_J3(t,XTheta,pConst), [t_im1 t_i], XTheta_full, opt);
    Theta_full = reshape(XTheta_full(end,n+1:end), n, 1);

        % Integrate Xstar, Phi, and Theta from t_im1 to t_i
    Phi_im1 = eye(n);
    XPhi_im1 = [Xstar_im1; reshape(Phi_im1,n^2,1)];
    [~, XPhi_i] = ode45(@(t,XPhi)STMEOM_J2(t,XPhi,pConst.mu, pConst.J2, pConst.Ri), [t_im1 t_i], XPhi_im1, opt);
    Xstar_i = XPhi_i(end,1:n)';
    Phi_i = reshape(XPhi_i(end,n+1:end),size(Phi_im1));

    Theta_im1 = zeros(n,1);
    XTheta_im1 = [Xstar_im1; reshape(Theta_im1,n,1)];
    [~, XTheta_i] = ode45(@(t,XTheta)STMEOM_CP_J3(t,XTheta,pConst), [t_im1 t_i], XTheta_im1, opt);
    Theta_i = reshape(XTheta_i(end,n+1:end),size(Theta_im1));

        % Build Psi(t, t_0)
    Psi_i = [
                Phi_full, Theta_full;
                zeros(1,n), eye(1)
            ];

        % Time update
    x_i = Phi_i*x_im1;
    P_i = Phi_i*P_im1*Phi_i';
    S_i = Phi_i*S_im1 + Theta_i;
    P_ci = P_i + S_i*Pcc0*S_i'; % aka Pxx^-
    P_xci = S_i*Pcc0;

        % Get number of measurements in Y, station states, and station 
        % visibility at this time
    meas = length(Y_i)/2; % Assuming 2 data points per measurement: range and range-rate
    Xstat = Xs{k}'; % Extract station state(s) at the time of measurement
    statVis = vis{k}; % Extract the stations that were visible at the time of measurement

        % Build y_i
    yExp = [];
    for kk = 1:meas
        genMeas = generateRngRngRate(Xstar_i, Xstat(:,meas), stations(statVis(kk)).elMask, true); % Ignore elevation mask
        yExp = [yExp; genMeas(1:2)];
    end

    y_i = Y_i - yExp;

        % Build Htilde_i
    Htilde_x_i = [];
    Htilde_c_i = [];
    for kk = 1:meas
        [Htilde_x, Htilde_c] = MeasurementPartials_CP_J3_sc(Xstar_i, Xstat(:,meas));
        Htilde_x_i = [Htilde_x_i; Htilde_x];
        Htilde_c_i = [Htilde_c_i; Htilde_c];
    end

        % Build K_i
    K_i = P_i*Htilde_x_i'*(Htilde_x_i*P_i*Htilde_x_i' + R_i)^-1;
    K_ci = P_ci*Htilde_x_i'*(Htilde_x_i*P_ci*Htilde_x_i' + R_i)^-1;

        % Measurement update
    xBar_i = x_i;
    x_i = x_i + K_i*(y_i - Htilde_x_i*x_i);
    S_i = (eye(n) - K_i*Htilde_x_i)*S_i - K_i*Htilde_c_i;
    x_ci = xBar_i + K_ci*(y_i - Htilde_x_i*x_i);

    mat = K_i*Htilde_x_i; % Intermediate matrix for sizing
    P_i = (eye(size(mat)) - mat)*P_i*(eye(size(mat)) - mat)' + K_i*R_i*K_i';
    P_ci = P_i + S_i*Pcc0'*S_i'; %% AKA Pxx^+
    P_xci = S_i*Pcc0;

        % Accumulate data to save
    xEst = [xEst, x_i];
    xcEst = [xcEst, x_ci];
    PEst = [PEst, {P_i}];
    PcEst = [PcEst, {P_ci}];
    PxcEst = [PxcEst, {P_xci}];
    prefit_res = [prefit_res, y_i];
    postfit_res = [postfit_res, y_i - Htilde_x_i*x_i];
    postfit_res_c = [postfit_res_c, y_i - Htilde_x_i*x_ci];
    XEst = [XEst, Xstar_i + x_i];
    XcEst = [XcEst, Xstar_i + x_ci];
    Phi_total = [Phi_total, {Phi_full}];
    Theta_total = [Theta_total, {Theta_full}];
    Psi = [Psi, {Psi_i}];

        % Update for next run
    t_im1 = t_i;
    Xstar_im1 = Xstar_i;
    P_im1 = P_i;
    x_im1 = x_i;
    S_im1 = S_i;

end

    %% Assign outputs
filterOut.xEst = xEst;
filterOut.xcEst = xcEst;
filterOut.PEst = PEst;
filterOut.PcEst = PcEst;
filterOut.PxcEst = PxcEst;
filterOut.prefit_res = prefit_res;
filterOut.postfit_res = postfit_res;
filterOut.postfit_res_c = postfit_res_c;
filterOut.t = t(2:end); % t_0 not included in estimate
filterOut.statVis = vis;
filterOut.XEst = XEst;
filterOut.XcEst = XcEst;
filterOut.Phi_total = Phi_total;
filterOut.Phi_test = Phi_test;
filterOut.Theta_total = Theta_total;
filterOut.Theta_test = Theta_test;
filterOut.Psi = Psi;

end