function filterOut = Potter(Xstar0, stations, pConst, scConst, P0, x0, numMeas)
% Function that implements the Potter algorithm for stat OD problems
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
%                   [xEst_0, xEst_1, ..., xEst_t-1], where 
%                   xEst = [x; y; z; xDot; yDot; zDot]
%           - PEst: Estimated state covariance at each time in t, organized 
%                   as follows:
%                   [{P_0}, {P_1}, ..., {P_t-1}]
%           - prefit_res: Pre-fit residuals (y_k) at each time in t:
%                         [y_0, y_1, ..., y_t-1]
%           - postfit_res: Postfit residuals (epsilon = y_k - Htilde_k*x_k)
%                          at each time in t:
%                          [epsilon_0, epsilon_1, ..., epsilon_t-1]
%           - t: Measurement time vector for the LKF filter
%           - statVis: Station visibility vector
%           - XEst: Estimated full state at each time in t:
%                   [XEst_0, XEst_1, ..., XEst_t-1], where
%                   XEst = [X; Y; Z; XDot; YDot; ZDot]
%           - Phi: Cell array of STMs from t0 to each t_k in t: 
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
Phi = eye(n);

    %% Process station data into a usable form
[t, Y, R, vis] = processStations(stations);

    % Specify number of measurements to process
if ~exist("numMeas", 'var')
    numMeas = length(Y);
end

    %% Loop through each observation
Xstar_k = Xstar0;
x_k = x0;
W_k = chol(P0); % Initialize W with cholesky decomposition of P0
Phi_full = eye(n);
for k = 1:numMeas
        % Read measurements and measurement uncertainty at this time
    t_k = t(k);
    Y_k = Y{k};
    R_k = R{k};
    statVis = vis{k};

        % Get number of measurements in Y
    meas = length(Y_k)/2; % Assuming 2 data points per measurement: range and range-rate

    res = [];
    for idx = 1:2
            % Choose measurement to process
        if idx == 1 % Process range
            measInclude = [true, false];
        else % Process range rate
            measInclude = [false, true];
        end

            % Build y_k
        yExp = [];
        for kk = 1:meas
            genMeas = generateRngRngRate(Xstar_k, statVis, stations(statVis(kk)).elMask, pConst, true); % Ignore elevation mask
            yExp = [yExp; genMeas(idx)];
        end
    
        y_k = Y_k(idx) - yExp;
        res = [res; y_k];

            % Build Htilde_k
        Htildes = [];
        Htilde_k = [];
        for kk = 1:meas
            Htilde_k = [Htilde_k; MeasurementPartials_RngRngRate(Xstar_k, statVis, pConst, measInclude)];
        end
        Htildes = [Htildes; Htilde_k];

            % Build Potter variables
        sigSqrd = R_k(idx,idx);
        Ftilde_k = W_k'*Htilde_k';
        alpha_k = (Ftilde_k'*Ftilde_k + sigSqrd)^-1;
        gamma_k = 1/(1+sqrt(alpha_k*sigSqrd));
        K_k = alpha_k*W_k*Ftilde_k;

            % Measurement update
        x_k = x_k + K_k*(y_k - Htilde_k*x_k);
        W_k = W_k - gamma_k*K_k*Ftilde_k';

    end

        % Accumulate data to save
    xEst = [xEst, x_k];
    PEst = [PEst, {W_k*W_k'}];
    prefit_res = [prefit_res, res];
    postfit_res = [postfit_res, res - Htildes*x_k];
    XEst = [XEst, Xstar_k + x_k];

    if k < numMeas
            % Read next time
        t_kp1 = t(k+1);
    
            % Continue to integrate Phi(t0, tf) for iteration purposes
        XPhi_full = [Xstar_k; reshape(Phi_full,n^2,1)];
        [~, XPhi_full] = ode45(@(t,XPhi)STMEOM_MuJ2Drag(t,XPhi,pConst,scConst), [t_k t_kp1], XPhi_full, opt);
        Phi_full = reshape(XPhi_full(end,n+1:end), n, n);
    
        Phi = [Phi; {Phi_full}];
    
            % Integrate Xstar and Phi from t_k to t_kp1
        Phi_k = eye(n);
        XPhi_k = [Xstar_k; reshape(Phi_k,n^2,1)];
        [~, XPhi_kp1] = ode45(@(t,XPhi)STMEOM_MuJ2Drag(t,XPhi,pConst,scConst), [t_k t_kp1], XPhi_k, opt);
        Xstar_kp1 = XPhi_kp1(end,1:n)';
        Phi_kp1 = reshape(XPhi_kp1(end,n+1:end),size(Phi_k));
    
            % Time update
        x_kp1 = Phi_kp1*x_k;
        W_kp1 = Phi_kp1*W_k;
    
            % Update for next run
        Xstar_k = Xstar_kp1;
        W_k = W_kp1;
        x_k = x_kp1;
    end

end



    %% Assign outputs
filterOut.xEst = xEst;
filterOut.PEst = PEst;
filterOut.prefit_res = prefit_res;
filterOut.postfit_res = postfit_res;
filterOut.t = t;
filterOut.statVis = vis;
filterOut.XEst = XEst;
filterOut.Phi = Phi;

end