function filterOut = SRIF(Xstar0, stations, pConst, P0, x0, Q0, uBar, forceUpperTriangular)
% Function that implements an SRIF for stat OD problems
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
%       - x0: Initial state deviation estimate
%       - Q0: Initial process noise covariance matrix
%       - uBar: Mean noise vector, generally zeros(3,1) but not always!
%       - forceUpperTriangular: Boolean indicating whether the time updated
%                               R matrix is forced to be upper triangular
%                               or not. Nominally, this should be true.
%   Outputs:
%       - filterOut: Filter output structure with the following fields:
%           - xEst: Estimated state deviation at each time processed from 
%                   the station measurements in "stations" (t), organized
%                   as follows:
%                   [xEst_1, xEst_2, ..., xEst_t], where 
%                   xEst = [x; y; z; xDot; yDot; zDot]
%           - PEst: Estimated state covariance at each time in t, organized 
%                   as follows:
%                   [{P_1}, {P_2}, ..., {P_t}]
%           - PBarEst: Estimated time update state covariance at each time
%                      in t, organized as follows:
%                      [{PBar_1}, {PBar_2}, ..., {PBar_t}];
%           - Ru: Square root noise covariance matrix used for smoothing:
%                 [{Ru_1}, {Ru_2}, ... {Ru_t}]
%           - Rux: Square root state noise covariance matrix for smoothing:
%                  [{Rux_1}, {Rux_2}, ..., {Rux_t}]
%           - bTildeu: Noise information vector for smoothing:
%                       [bTildeu_1, bTildeu_2, ..., bTildeu_t]
%           - uHat: Process noise vector at each time in t until t_tm1:
%                   [uHat_0, uHat_1, ..., uHat_tm1]
%           - prefit_res_whitened: Whitened pre-fit residuals (y_i) at each
%                                  time in t:
%                                  [y_1, y_2, ..., y_t]
%           - postfit_res_whitened: Whitened postfit residuals at each time
%                                   in t:
%                                   [e_1, e_2, ..., e_t]
%           - prefit_res: Un-whitened pre-fit residuals (V_i*y_i) at each
%                         time in t:
%                         [V_1*y_1, V_2*y_2, ..., V_t*y_t]
%           - postfit_res: Un-whitened postfit residuals at each time in t:
%                          [V_1*e_1, V_2*e_2, ..., V_t*e_t]
%           - t: Measurement time vector for the LKF filter
%           - statVis: Station visibility vector
%           - XStar: Nominal full state at each time in t:
%                    [XStar_1, XStar_2, ..., XStar_t]
%           - XEst: Estimated full state at each time in t, computed as 
%                   XEst = XNom + xEst:
%                   [XEst_1, XEst_2, ..., XEst_t], where
%                   XEst = [X; Y; Z; XDot; YDot; ZDot]
%           - Phi_total: Cell array of STMs from t0 to each t_i in t: 
%                        [{Phi(t_1, t0)}; {Phi(t_2,t0)}; ...; {Phi(t_f,t_0)}]
%           - Phi: Cell array of STMs from t_im1 to t_i:
%                  [{Phi(t_1, t_0)}; {Phi(t_2,t_1)}; ...; {Phi(t_i,t_im1}]
%
%   By: Ian Faber, 03/07/2025
%

%% Initialize settings

    % Default to forcing an upper triangular R
if isempty(forceUpperTriangular)
    forceUpperTriangular = true;
end

    % Format ode45
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);

    % Find state and noise sizes
n = length(Xstar0); % Length of state
q = length(uBar); % Length of noise

    % Preallocate filter outputs
xEst = [];
PEst = [];
PBarEst = [];
Ru = [];
Rux = [];
bTildeu = [];
uHat = [];
prefit_res_whitened = [];
postfit_res_whitened = [];
prefit_res = [];
postfit_res = [];
XStar = [];
XEst = [];
Phi_total = [];
Phi = [];

%% Define helper function
GammaFunc = @(dt) [(dt/2)*eye(3); eye(3)];

%% Process station data into a usable form
[t, Y, measCov, Xs, vis] = processStations(stations); % measCov = R in a normal Kalman Filter!

%% Whiten observations
for k = 1:length(Y)
    V = chol(measCov{k}, 'lower'); % measCov = V*V'
    
    Y{k} = (V^-1)*Y{k};
end

%% Loop through each observation
t_im1 = t(1);
Xstar_im1 = Xstar0;
x_im1 = x0;
R_im1 = chol(P0^-1,"upper"); % Find Rbar_0 -> P0^-1 = Lambda = R'*R
if any(Q0 > 0)
    Ru_im1 = chol(Q0^-1,"upper"); % Find Ru_0 -> Q = (Ru'*Ru)^-1
    bu_im1 = Ru_im1*uBar;
end
b_im1 = R_im1*x_im1; % Find bbar_0
Phi_full = eye(n);
for k = 2:length(Y)
        % Read next time, measurement, and measurement covariance
    t_i = t(k);
    Y_i = Y{k};
    V_i = chol(measCov{k},'lower');

        % Continue to integrate Phi(t0, tf) for iteration purposes
    XPhi_full = [Xstar_im1; reshape(Phi_full,n^2,1)];
    [~, XPhi_full] = ode45(@(t,XPhi)STMEOM_J2(t,XPhi,pConst.mu, pConst.J2, pConst.Ri), [t_im1 t_i], XPhi_full, opt);
    Phi_full = reshape(XPhi_full(end,n+1:end), n, n);

    Phi_total = [Phi_total; {Phi_full}];

        % Integrate Xstar and Phi from t_im1 to t_i
    Phi_im1 = eye(n);
    XPhi_im1 = [Xstar_im1; reshape(Phi_im1,n^2,1)];
    [~, XPhi_i] = ode45(@(t,XPhi)STMEOM_J2(t,XPhi,pConst.mu, pConst.J2, pConst.Ri), [t_im1 t_i], XPhi_im1, opt);
    Xstar_i = XPhi_i(end,1:n)';
    Phi_i = reshape(XPhi_i(end,n+1:end),size(Phi_im1)); % Phi(t_i, t_im1)

    Phi = [Phi; {Phi_i}];

        % Time update
    delT = t_i - t_im1;
    if any(any(Q0 > 0) & (delT <= 10)) % Time update with process noise
            % Make Gamma and Q
        Gamma_i = GammaFunc(delT);
        Q_i = Q0;
        
            % Set up noise variables
        Ru_k = chol(Q_i^-1,"upper"); % Q = (Ru'*Ru)^-1        

            % Set up Rtilde
        Rtilde = R_im1*(Phi_i^-1);
        
        mat = [
                Ru_k, zeros(q,n), bu_im1;
                -Rtilde*Gamma_i, Rtilde, b_im1
              ];

            % Run Householder
        out = Householder(mat);

            % Extract time update outputs
        R_i = out(q+1:end, q+1:q+n);
        b_i = out(q+1:end, q+n+1);

            % Extract process noise outputs
        Ru_k = out(1:q,1:q);
        Rux_k = out(1:q,q+1:q+n);
        bTildeu_k = out(1:q,q+n+1);

            % Propagate noise information vector
        bu_im1 = Ru_k*uBar;

            % Accumulate process noise outputs
        Ru = [Ru, {Ru_k}];
        Rux = [Rux, {Rux_k}];
        bTildeu = [bTildeu, bTildeu_k];

    else  % Time update without process noise
        x_i = Phi_i*x_im1;
        R_i = R_im1*(Phi_i^-1);
        % b_i = R_i*x_i; % This should be constant across a time update! b_i = R_i*x_i = R_im1*(Phi_i^-1)*Phi_i*x_im1 = R_im1*x_im1 = b_im1
        b_i = b_im1;

            % Force R_i to be upper triangular with Householder Transformation
        if forceUpperTriangular
            out = Householder([R_i,b_i]);
            R_i = out(1:n,1:n);
            b_i = out(1:n,n+1);
        end
    end
    
        % Reconstruct PBar
    PBarEst = [PBarEst, {(R_i'*R_i)^-1}]; % P = (Lambda)^-1 = (R'*R)^-1

        % Get number of measurements in Y, station states, and station 
        % visibility at this time
    meas = length(Y_i)/2; % Assuming 2 data points per measurement: range and range-rate
    Xstat = Xs{k}'; % Extract station state(s) at the time of measurement
    statVis = vis{k}; % Extract the stations that were visible at the time of measurement

        % Build y_i
    yExp = [];
    for kk = 1:meas
        genMeas = generateRngRngRate(Xstar_i, Xstat(:,meas), stations(statVis(kk)).elMask, true); % Ignore elevation mask
        yExp = [yExp; (V_i^-1)*genMeas(1:2)]; % Whiten the expected measurement
    end

    y_i = Y_i - yExp;

        % Build Htilde_i
    Htilde_i = [];
    for kk = 1:meas
        Htilde = MeasurementPartials_RngRngRate_sc(Xstar_i, Xstat(:,meas));
        Htilde_i = [Htilde_i; (V_i^-1)*Htilde]; % Whiten the Htilde matrix
    end

        % Measurement update
    mat = [R_i, b_i; Htilde_i, y_i];
    out = Householder(mat);

    R_i = out(1:n,1:n);
    b_i = out(1:n,n+1);
    e = out(n+1:end,n+1);

        % Intermediate variable for xHat
    xHat = (R_i^-1)*b_i;

    if any(any(Q0>0) & (delT <= 10))
            % Calculate uHat_im1 and accumulate it
        u_im1 = (Ru_k^-1)*(bTildeu_k - Rux_k*xHat);
        uHat = [uHat, u_im1];
    end

        % Accumulate data to save
    xEst = [xEst, xHat]; % xHat = R^-1*b
    PEst = [PEst, {(R_i'*R_i)^-1}]; % P = (Lambda)^-1 = (R'*R)^-1
    prefit_res_whitened = [prefit_res_whitened, y_i];
    postfit_res_whitened = [postfit_res_whitened, e];
    prefit_res = [prefit_res, V_i*y_i];
    postfit_res = [postfit_res, V_i*e];
    XStar = [XStar, Xstar_i];
    XEst = [XEst, Xstar_i + xHat];

        % Update for next run
    t_im1 = t_i;
    Xstar_im1 = Xstar_i;
    x_im1 = xHat;
    R_im1 = R_i;
    b_im1 = b_i;

end

%% Assign outputs
filterOut.xEst = xEst;
filterOut.PEst = PEst;
filterOut.PBarEst = PBarEst;
filterOut.Ru = Ru;
filterOut.Rux = Rux;
filterOut.bTildeu = bTildeu;
filterOut.uHat = uHat;
filterOut.prefit_res_whitened = prefit_res_whitened;
filterOut.postfit_res_whitened = postfit_res_whitened;
filterOut.prefit_res = prefit_res;
filterOut.postfit_res = postfit_res;
filterOut.t = t(2:end); % t_0 not included in estimate
filterOut.statVis = vis;
filterOut.XStar = XStar;
filterOut.XEst = XEst;
filterOut.Phi_total = Phi_total;
filterOut.Phi = Phi;


end