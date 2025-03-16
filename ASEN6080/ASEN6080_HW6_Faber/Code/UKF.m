function filterOut = UKF(stations, pConst, X0, P0, Q0, alpha, beta, includeJ3)
% Function that implements a UKF for stat OD problems
%   Inputs:
%       - stations: Stations struct as defined by makeStations.m. Must have
%                   propagated station states! To propagate states, see
%                   generateTruthData.m.
%       - pConst: Planetary constant structure as formatted by
%                 getPlanetConst.m
%       - X0: Initial full state estimate
%       - P0: Initial state covariance estimate
%       - Q0: Initial process noise covariance matrix
%       - alpha: UKF sigma point spacing variable from [1e-4, 1]
%       - beta: UKF probability distribution variable, generally 2 for 
%               Gaussian probability distributions
%       - includeJ3: Boolean indicating whether the filter dynamics should
%                    include J3 in addition to mu and J2
%   Outputs:
%       - filterOut: Output filter structure with the following fields:
%           - XEst: Estimated full state vector at each time in t:
%                   [XEst_1, XEst_2, ..., XEst_t], where
%                   XEst = [X; Y; Z; XDot; YDot; ZDot]
%           - PEst: Estimated state covariance at each time in t, organized 
%                   as follows:
%                   [{P_1}, {P_2}, ..., {P_t}]
%           - prefit_res: Pre-fit residuals (y_i - yBar_i) at each time in t:
%                         [prefit_1, prefit_2, ..., prefit_t]
%           - postfit_res: Post-fit residuals (y_i - yBar_i after XEst has 
%                          been computed) at each time in t:
%                          [postfit_1, postfit_2, ..., postfit_t]
%           - t: Measurement time vector for the EKF filter
%           - statVis: Station visibility vector
%
%   By: Ian Faber, 03/15/2025
%
%% Initialize settings
    % Format ode45 and sizes
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
L = length(X0);
XEst = [];
PEst = [];
prefit_res = [];
postfit_res = [];

%% Define helper functions
GammaFunc = @(dt) [(dt/2)*eye(3); eye(3)];
J2Func = @(t,X)orbitEOM_MuJ2(t,X,pConst.mu,pConst.J2,pConst.Ri);
J3Func = @(t,X)orbitEOM_MuJ2J3(t,X,pConst.mu,pConst.J2,pConst.J3,pConst.Ri);

%% Process station data into a usable form
[t, Y, R, Xs, vis] = processStations(stations);

%% Precompute UKF weights
kappa = 3 - L;
lambda = (alpha^2)*(L + kappa) - L;
gamma = sqrt(L + lambda);

W_0m = lambda/(L + lambda);
W_0c = lambda/(L + lambda) + (1 - alpha^2 + beta);
W_im = [W_0m, (1/(2*(L + lambda)))*ones(1,2*L)];
W_ic = [W_0c, (1/(2*(L + lambda)))*ones(1,2*L)];

%% Loop through all observations
    % Initialize UKF variables
X_im1 = X0;
P_im1 = P0;
t_im1 = t(1);

    % Loop through each observation
for k = 2:length(Y)
        % Read next time, measurement, and measurement covariance
    t_i = t(k);
    Y_i = Y{k};
    R_i = R{k};

        % Create Q if necessary for process noise
    dT = t_i - t_im1;
    if any(any(Q0 > 0) & (dT <= 10)) % Process noise exists and the time gap isn't too big for assumptions to break
        Gamma_i = GammaFunc(dT);
        Q = Gamma_i*Q0*Gamma_i';
    else
        Q = zeros(L);
    end

        % Calculate previous sigma points
    sqrtP_im1 = sqrtm(P_im1); % Used to be chol()
    Chi_im1 = [X_im1, X_im1 + gamma*sqrtP_im1, X_im1 - gamma*sqrtP_im1]; % L x (2L + 1) matrix

        % Propagate previous sigma points through dynamics
    ChiVec_im1 = reshape(Chi_im1, L*(2*L+1), 1);
    tspan = [t_im1, t_i];
    if ~includeJ3 % Only include mu and J2
        [~,ChiVec] = ode45(@(t,ChiVec)sigPointEOM(t,ChiVec,J2Func), tspan, ChiVec_im1, opt);
    else % Include mu, J2, and J3
        [~,ChiVec] = ode45(@(t,ChiVec)sigPointEOM(t,ChiVec,J3Func), tspan, ChiVec_im1, opt);
    end
    Chi_i = reshape(ChiVec(end,:), L, 2*L + 1);

        % Time update
    X_i = 0;
    for kk = 1:2*L+1
        X_i = X_i + W_im(kk)*Chi_i(:,kk);
    end
    
    P_i = Q;
    for kk = 1:2*L+1
        P_i = P_i + W_ic(kk)*(Chi_i(:,kk) - X_i)*(Chi_i(:,kk) - X_i)';
    end

        % Recompute sigma points to account for propagation and process
        % noise
    sqrtP_i = sqrtm(P_i); % Used to be chol()
    Chi_i = [X_i, X_i + gamma*sqrtP_i, X_i - gamma*sqrtP_i];

        % Get number of measurements in Y, station states, and station 
        % visibility at this time
    meas = length(Y_i)/2; % Assuming 2 data points per measurement: range and range-rate
    Xstat = Xs{k}'; % Extract station state(s) at the time of measurement
    statVis = vis{k}; % Extract the stations that were visible at the time of measurement
    
        % Construct yBar_i
    yBar_i = 0;
    YExp = [];
    for kk = 1:2*L + 1
        yExp = [];
        state = Chi_i(:,kk);
        for idx = 1:meas % Account for multiple stations visible at the same time
            genMeas = generateRngRngRate(state, Xstat(:,idx), stations(statVis(idx)).elMask, true); % Ignore elevation mask
            yExp = genMeas(1:2);
            % YExp = [YExp, genMeas];
            YExp = [YExp, yExp];
            yBar_i = yBar_i + W_im(kk)*yExp;
        end
    end
    
        % Compute innovation and cross covariances
    Pyy = R_i;
    Pxy = zeros(L,2);
    for kk = 1:2*L + 1
        Pyy = Pyy + W_ic(kk)*(YExp(:,kk) - yBar_i)*(YExp(:,kk) - yBar_i)';
        Pxy = Pxy + W_ic(kk)*(Chi_i(:,kk) - X_i)*(YExp(:,kk) - yBar_i)';
    end

        % Compute Kalman Gain
    K_i = Pxy*(Pyy^-1);

        % Measurement update
    X_i = X_i + K_i*(Y_i - yBar_i);
    P_i = P_i - K_i*Pyy*K_i';
    
        % Calculate expected measurement after measurement update for
        % postfits
    genMeas_post = generateRngRngRate(X_i, Xstat(:,1), stations(statVis(1)).elMask, true);

        % Accumulate data to save
    XEst = [XEst, X_i];
    PEst = [PEst, {P_i}];
    prefit_res = [prefit_res, Y_i - yBar_i];
    postfit_res = [postfit_res, Y_i - genMeas_post(1:2)];

        % Update for next run
    X_im1 = X_i;
    P_im1 = P_i;
    t_im1 = t_i;

end

%% Assign outputs
filterOut.XEst = XEst;
filterOut.PEst = PEst;
filterOut.prefit_res = prefit_res;
filterOut.postfit_res = postfit_res;
filterOut.t = t(2:end); % t_0 not included in estimate
filterOut.statVis = vis;


end