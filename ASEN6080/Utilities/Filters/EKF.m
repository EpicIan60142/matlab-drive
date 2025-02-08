function filterOut = EKF(Xstar0, stations, pConst, P0, t_start)
% Function that implements an EKF for stat OD problems
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
%       - t_start: Time to start the filter, generally the time after some
%                  number of LKF measurements have been taken
%   Outputs:
%       - filterOut: Filter output structure with the following fields:
%           - xEst: Estimated state deviation at each time processed from 
%                   the station measurements in "stations" (t), organized
%                   as follows:
%                   [xEst_1, xEst_2, ..., x_t], where 
%                   xEst = [x; y; z; xDot; yDot; zDot]
%           - PEst: Estimated state covariance at each time in t, organized 
%                   as follows:
%                   [{P_1}, {P_2}, ..., {P_t}]
%           - postfit_res: Post-fit residuals (y_i) at each time in t:
%                          [y_1, y_2, ..., y_t]
%           - epsilon: Observation errors (y_i - Htilde_i*x_i) at each time
%                      in t:
%                      [epsilon_1, epsilon_2, ..., epsilon_t]
%           - t: Measurement time vector for the EKF filter
%           - statVis: Station visibility vector
%           - XEst: Estimated full state at each time in t:
%                   [XEst_1, XEst_2, ..., XEst_t], where
%                   XEst = [X; Y; Z; XDot; YDot; ZDot]
%           - Phi_full: Integrated STM from t0 to tf, for iteration
%                       purposes (Phi(t0, tf))
%
%   By: Ian Faber, 02/02/2025
%

%% Initialize settings

    % Format ode45 and sizes
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
n = length(Xstar0);
xEst = [];
PEst = [];
postfit_res = [];
epsilon = [];
XEst = [];

    %% Process station data into a usable form
[t, Y, R, Xs, vis] = processStations(stations, t_start);

    %% Loop through each observation
t_im1 = t(1);
Xstar_im1 = Xstar0;
P_im1 = P0;
for k = 2:length(Y)
    
        % Read next time, measurement, and measurement covariance
    t_i = t(k);
    Y_i = Y{k};
    R_i = R{k};

        % Integrate Xstar and Phi from t_im1 to t_i
    Phi_im1 = eye(n);
    XPhi_im1 = [Xstar_im1; reshape(Phi_im1,n^2,1)];
    [~, XPhi_i] = ode45(@(t,XPhi)STMEOM_J2(t,XPhi,pConst.mu, pConst.J2, pConst.Ri), [t_im1 t_i], XPhi_im1, opt);
    Xstar_i = XPhi_i(end,1:n)';
    Phi_i = reshape(XPhi_i(end,n+1:end),size(Phi_im1));

        % Time update
    P_i = Phi_i*P_im1*Phi_i';

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
    Htilde_i = [];
    for kk = 1:meas
        Htilde_i = [Htilde_i; MeasurementPartials_RngRngRate_sc(Xstar_i, Xstat(:,meas))];
    end

        % Build K_i
    K_i = P_i*Htilde_i'*(Htilde_i*P_i*Htilde_i' + R_i)^-1;
    
        % Measurement and reference orbit update
    x_i = K_i*y_i;
    Xstar_i = Xstar_i + x_i;

    mat = K_i*Htilde_i; % Intermediate matrix for sizing
    P_i = (eye(size(mat)) - mat)*P_i*(eye(size(mat)) - mat)' + K_i*R_i*K_i';

        % Accumulate data to save
    xEst = [xEst, x_i];
    PEst = [PEst, {P_i}];
    postfit_res = [postfit_res, y_i];
    epsilon = [epsilon, y_i - Htilde_i*x_i];
    XEst = [XEst, Xstar_i];

        % Update for next run
    t_im1 = t_i;
    Xstar_im1 = Xstar_i;
    P_im1 = P_i;

end

    %% Assign outputs
filterOut.xEst = xEst;
filterOut.PEst = PEst;
filterOut.postfit_res = postfit_res;
filterOut.epsilon = epsilon;
filterOut.t = t(2:end); % t_0 not included in estimate
filterOut.statVis = vis;
filterOut.XEst = XEst;

end