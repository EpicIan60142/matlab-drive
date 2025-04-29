function filterOut = IEKF_DMC(Xstar0, stations, pConst, scConst, P0, B, Qu, t_start, t_end)
% Function that implements an iterated EKF for stat OD problems
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
%           - prefit_res: Pre-fit residuals (y_i) at each time in t:
%                         [y_1, y_2, ..., y_t]
%           - postfit_res: Post-fit residuals (epsilon = y_i - Htilde_i*x_i)
%                          at each time in t:
%                          [epsilon_1, epsilon_2, ..., epsilon_t]
%           - t: Measurement time vector for the EKF filter
%           - statVis: Station visibility vector
%           - XEst: Estimated full state at each time in t:
%                   [XEst_1, XEst_2, ..., XEst_t], where
%                   XEst = [X; Y; Z; XDot; YDot; ZDot]s
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

    %% Process station data into a usable form
[t, Y, R, Xs, vis] = processStations(stations, t_start);

if exist("t_end", 'var')
    kEnd = find(t == t_end, 1, 'first');
else
    kEnd = length(Y);
end
    %% Loop through each observation
t_im1 = t(1);
Xstar_im1 = Xstar0;
P_im1 = P0;
for k = 2:kEnd
    
        % Read next time, measurement, and measurement covariance
    t_i = t(k);
    Y_i = Y{k};
    R_i = R{k};

        % Generate white gaussian noise
    u = randn(3,1);
    u = chol(Qu)*u; % Scale noise properly

        % Integrate Xstar and Phi from t_im1 to t_i
    Phi_im1 = eye(n);
    XPhi_im1 = [Xstar_im1; reshape(Phi_im1,n^2,1)];
    [~, XPhi_i] = ode45(@(t,XPhi)STMEOM_MuSunSRP_DMC(t,XPhi,B,u,pConst,scConst), [t_im1 t_i], XPhi_im1, opt);
    Xstar_i = XPhi_i(end,1:n)';
    Phi_i = reshape(XPhi_i(end,n+1:end),size(Phi_im1));

        % Make Q_i
    Q_i = makeQ_DMC(B, Qu, t_i, t_im1);

        % Adjust for C_R
    Q_iX = Q_i(1:6,1:6);
    Q_iW_1 = Q_i(1:6,7:9);
    Q_iW_2 = Q_i(7:9,1:6);
    Q_iW_3 = Q_i(7:9,7:9);

    Q_i = [
            Q_iX, zeros(6,1), Q_iW_1;
            zeros(1,6), 0, zeros(1,3);
            Q_iW_2, zeros(3,1), Q_iW_3
          ];

    if t_i - t_im1 > 60 % Turn off DMC
        Q_i = 0*Q_i;
    end

        % Time update
    P_i = Phi_i*P_im1*Phi_i' + Q_i;

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
        Htilde_i = [Htilde_i; MeasurementPartials_RngRngRate_sc_DMC(Xstar_i, Xstat(:,meas))];
    end

        % Build K_i
    K_i = P_i*Htilde_i'*(Htilde_i*P_i*Htilde_i' + R_i)^-1;

        % Measurement and reference orbit update
    x_im = zeros(size(Xstar_i));
    x_i = K_i*y_i;
    Xstar_i = Xstar_i + x_i;

    mat = K_i*Htilde_i; % Intermediate matrix for sizing
    P_i = (eye(size(mat)) - mat)*P_i*(eye(size(mat)) - mat)' + K_i*R_i*K_i';

        % Repeat measurement update until residuals reach an acceptable
        % multiple of sigma
    iter = 0;
    boundLevel = 2;
    while ((abs(y_i(1)) > boundLevel*sqrt(R_i(1,1))) || (abs(y_i(2)) > boundLevel*sqrt(R_i(2,2))))
        if iter >= 9999
            % fprintf("\n\t\tIEKF iterated %.0f times..., y_i = [%.3f, %.3f]", iter, y_i)
            % fprintf("\nBreaking out")
            break;
        end
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
            Htilde_i = [Htilde_i; MeasurementPartials_RngRngRate_sc_DMC(Xstar_i, Xstat(:,meas))];
        end

            % Build K_i
        K_i = P_i*Htilde_i'*(Htilde_i*P_i*Htilde_i' + R_i)^-1;

            % Measurement and reference orbit update
        x_i = x_im + K_i*(y_i - Htilde_i*(x_im-x_i));
        Xstar_i = Xstar_i + x_i;

        mat = K_i*Htilde_i; % Intermediate matrix for sizing
        P_i = (eye(size(mat)) - mat)*P_i*(eye(size(mat)) - mat)' + K_i*R_i*K_i';

        iter = iter + 1;
    end

    if iter > 0
        if iter == 1
            fprintf("\n\t\tt = %.0f sec (%.3f days): IEKF iterated %.0f time! y_i = [%.3e, %.3e]", t_i, t_i/(24*60*60), iter, y_i)
        else
            fprintf("\n\t\tt = %.0f sec (%.3f days): IEKF iterated %.0f times! y_i = [%.3e, %.3e]", t_i, t_i/(24*60*60), iter, y_i)
        end
    end

        % Accumulate data to save
    xEst = [xEst, x_i];
    PEst = [PEst, {P_i}];
    prefit_res = [prefit_res, y_i];
    postfit_res = [postfit_res, y_i - Htilde_i*x_i];
    XEst = [XEst, Xstar_i];

        % Update for next run
    t_im1 = t_i;
    Xstar_im1 = Xstar_i;
    P_im1 = P_i;

end

    %% Assign outputs
filterOut.xEst = xEst;
filterOut.PEst = PEst;
filterOut.prefit_res = prefit_res;
filterOut.postfit_res = postfit_res;
filterOut.t = t(2:kEnd); % t_0 not included in estimate
filterOut.statVis = vis;
filterOut.XEst = XEst;

end