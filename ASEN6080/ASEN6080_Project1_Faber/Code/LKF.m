function filterOut = LKF(Xstar0, stations, pConst, scConst, P0, x0, numMeas)
% Function that implements an LKF for stat OD problems
%   Inputs:
%       - Xstar0: Initial value of the reference trajectory, organized as
%                 follows:
%                 [X0; Y0; Z0; Xdot0; Ydot0; Zdot0]
%       - stations: Stations struct as defined by makeStations.m. Must have
%                   propagated station states! To propagate states, see
%                   generateTruthData.m.
%       - pConst: Planetary constant structure as formatted by
%                 getPlanetConst.m
%       - scConst: Spacecraft constant structure as formatted by
%                  getSCConst.m
%       - P0: Initial state covariance estimate
%       - x0: Initial state deviation estimate
%       - numMeas: Number of measurements to process (optional). If not 
%                  specified, defaults to all measurements
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
%           - prefit_res: Pre-fit residuals (y_i) at each time in t:
%                         [y_1, y_2, ..., y_t]
%           - postfit_res: Postfit residuals (epsilon = y_i - Htilde_i*x_i)
%                          at each time in t:
%                          [epsilon_1, epsilon_2, ..., epsilon_t]
%           - t: Measurement time vector for the LKF filter
%           - statVis: Station visibility vector
%           - XEst: Estimated full state at each time in t:
%                   [XEst_1, XEst_2, ..., XEst_t], where
%                   XEst = [X; Y; Z; XDot; YDot; ZDot]
%           - Phi: Cell array of STMs from t0 to each t_i in t: 
%                  [{Phi(t_1, t0)}; {Phi(t_2,t0)}; ...; {Phi(t_f,t_0)}]
%
%   By: Ian Faber, 02/02/2025
%

%% Initialize settings

    % Format ode45 and sizes
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
n = length(Xstar0);
xEst = [];
PEst = [];
prefit_res = [];
postfit_res = [];
XEst = [];
Phi = [];

    %% Process station data into a usable form
[t, Y, R, vis] = processStations(stations);

    % Specify number of measurements to process
if ~exist("numMeas", 'var')
    numMeas = length(Y);
end

    %% Loop through each observation
t_im1 = t(1);
Xstar_im1 = Xstar0;
x_im1 = x0;
P_im1 = P0;
Phi_full = eye(n);
for k = 2:numMeas

        % Read next time, measurement, and measurement covariance
    t_i = t(k);
    Y_i = Y{k};
    R_i = R{k};

        % Continue to integrate Phi(t0, tf) for iteration purposes
    XPhi_full = [Xstar_im1; reshape(Phi_full,n^2,1)];
    [~, XPhi_full] = ode45(@(t,XPhi)STMEOM_MuJ2Drag(t,XPhi,pConst,scConst), [t_im1 t_i], XPhi_full, opt);
    Phi_full = reshape(XPhi_full(end,n+1:end), n, n);

    Phi = [Phi; {Phi_full}];

        % Integrate Xstar and Phi from t_im1 to t_i
    Phi_im1 = eye(n);
    XPhi_im1 = [Xstar_im1; reshape(Phi_im1,n^2,1)];
    [~, XPhi_i] = ode45(@(t,XPhi)STMEOM_MuJ2Drag(t,XPhi,pConst,scConst), [t_im1 t_i], XPhi_im1, opt);
    Xstar_i = XPhi_i(end,1:n)';
    Phi_i = reshape(XPhi_i(end,n+1:end),size(Phi_im1));

        % Time update
    x_i = Phi_i*x_im1;
    P_i = Phi_i*P_im1*Phi_i';
    % chol(P_i)

        % Get number of measurements in Y, station states, and station 
        % visibility at this time
    meas = length(Y_i)/2; % Assuming 2 data points per measurement: range and range-rate
    statVis = vis{k}; % Extract the stations that were visible at the time of measurement

        % Build y_i
    yExp = [];
    for kk = 1:meas
        genMeas = generateRngRngRate(Xstar_i, statVis, stations(statVis(kk)).elMask, pConst, true); % Ignore elevation mask
        yExp = [yExp; genMeas(1:2)];
    end

    y_i = Y_i - yExp;

        % Build Htilde_i
    Htilde_i = [];
    for kk = 1:meas
        Htilde_i = [Htilde_i; MeasurementPartials_RngRngRate(Xstar_i, statVis, pConst)];
    end

        % Build K_i
    K_i = P_i*Htilde_i'*(Htilde_i*P_i*Htilde_i' + R_i)^-1;

        % Measurement update
    x_i = x_i + K_i*(y_i - Htilde_i*x_i);

    mat = K_i*Htilde_i; % Intermediate matrix for sizing
    P_i = (eye(size(mat)) - mat)*P_i*(eye(size(mat)) - mat)' + K_i*R_i*K_i';
    % P_i = (eye(size(mat))-mat)*P_i;

        % Accumulate data to save
    xEst = [xEst, x_i];
    PEst = [PEst, {P_i}];
    prefit_res = [prefit_res, y_i];
    postfit_res = [postfit_res, y_i - Htilde_i*x_i];
    XEst = [XEst, Xstar_i + x_i];

        % Update for next run
    t_im1 = t_i;
    Xstar_im1 = Xstar_i;
    P_im1 = P_i;
    x_im1 = x_i;

end

    %% Assign outputs
filterOut.xEst = xEst;
filterOut.PEst = PEst;
filterOut.prefit_res = prefit_res;
filterOut.postfit_res = postfit_res;
filterOut.t = t(2:end); % t_0 not included in estimate
filterOut.statVis = vis;
filterOut.XEst = XEst;
filterOut.Phi = Phi;

end