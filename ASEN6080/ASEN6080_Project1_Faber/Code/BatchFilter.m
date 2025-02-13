function filterOut = BatchFilter(Xstar0, stations, pConst, scConst, P0, x0, dt)
% Function that implements a Batch Filter for OD problems.
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
%       - P0: Initial state covariance estimate (optional)
%             - Ignore by passing in [] instead
%       - x0: Initial state deviation estimate (optional)
%             - Ignore by passing in [] instead
%       - dt: Timestep for outputing integrated STMs (optional)
%             - Ignore by passing in [] instead
%   Outputs:
%       - filterOut: Filter output structure with the following fields:
%           - x0Est: Estimated initial state deviation, organized as follows:
%                   [x0; y0; z0; xDot0; yDot0; zDot0]
%           - P0Est: Estimated initial state covariance
%           - prefit_res: Pre-fit residuals (y_i) at each time in t:
%                         [y_1, y_2, ..., y_t]
%           - postfit_res: Post-fit residuals (epsilon_i = y_i - H_i*x_0)
%                          at each time in t:
%                          [epsilon_1, epsilon_2, ..., epsilon_t]
%           - t: Measurement time vector for the batch filter
%           - statVis: Station visibility vector
%           - Phi: Cell array of STMs from t0 to each t_i in t: 
%                  [{Phi(t_1, t0)}; {Phi(t_2,t0)}; ...; {Phi(t_f,t_0)}]
%
%   By: Ian Faber, 01/30/2025
%
    %% Initialize variables based on what's been provided
if ~isempty(P0) && ~isempty(x0)
    Lambda = P0^-1;
    N = Lambda*x0;
else
    Lambda = zeros(length(Xstar0));
    N = zeros(size(Xstar0));
end

timeStep = false; % Boolean indicating if a constant output timestep exists
if ~isempty(dt)
    timeStep = true;
end

    % Format ode45 and sizes
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
n = length(Xstar0);
prefit_res = [];
postfit_res = [];
Phi = [];

    %% Process station data into a usable form
[t, Y, R, vis] = processStations(stations);
Hs = [];

    %% Loop through each observation
t_im1 = t(1); % Start t at t_0
Xstar_im1 = Xstar0; % Start Xstar at Xstar(t_0)
Phi_im1 = eye(length(Xstar0));
for k = 2:length(Y)

        % Read next time, measurement, and measurement covariance
    t_i = t(k);
    Y_i = Y{k};
    R_i = R{k};

        % Integrate STM and EOM from t_{i-1} to t_i
    if timeStep % If a constant output timestep is defined, use it
        tspan = t_im1:dt:t_i;
    else
        tspan = [t_im1 t_i];
    end

    XPhi_im1 = [Xstar_im1; reshape(Phi_im1,n^2,1)];
    [t_int, XPhi_i] = ode45(@(t,XPhi)STMEOM_MuJ2Drag(t,XPhi,pConst,scConst), tspan, XPhi_im1, opt);
    Xstar_i = XPhi_i(end,1:n)';
    Phi_i = reshape(XPhi_i(end,n+1:end),size(Phi_im1));

    if timeStep
        for kk = 2:length(tspan)
            Phi = [Phi; {reshape(XPhi_i(t_int == tspan(kk), n+1:end), size(Phi_i)), tspan(kk)}];
            % kk
        end
    else
        Phi = [Phi; {Phi_i, tspan(2)}];
    end

        % Build Htilde_i
    meas = length(Y_i)/2; % Find number of measurements in this observation
    Htilde_i = [];
    statVis = vis{k}; % Extract the stations that were visible at the time of measurement
    for kk = 1:meas
        Htilde_i = [Htilde_i; MeasurementPartials_RngRngRate(Xstar_i, statVis, pConst)];
    end

        % Build y_i
    yExp = [];
    for kk = 1:meas
        genMeas = generateRngRngRate(Xstar_i, statVis, stations(statVis(kk)).elMask, pConst, true); % Ignore elevation mask
        yExp = [yExp; genMeas(1:2)]; % Only pull out range and range rate
    end

    y_i = Y_i - yExp;

    prefit_res = [prefit_res, y_i(1:2,:)];

        % Build H_i
    H_i = Htilde_i*Phi_i;
    Hs = [Hs; {H_i}];

        % Accumulate observation
    Lambda = Lambda + H_i'*(R_i^-1)*H_i;
    N = N + H_i'*(R_i^-1)*y_i;

        % Update for next iteration
    t_im1 = t_i;
    Xstar_im1 = Xstar_i;
    Phi_im1 = Phi_i;

end

    %% Solve normal equations
P0Est = Lambda^-1;
x0Est = P0Est*N;

    %% Assign epsilon
for k = 1:length(Y)-1
    postfit_res(:,k) = prefit_res(:,k) - Hs{k}*x0Est;
end

    %% Assign outputs
filterOut.x0Est = x0Est;
filterOut.P0Est = P0Est;
filterOut.prefit_res = prefit_res;
filterOut.postfit_res = postfit_res;
filterOut.t = t(2:end); % t_0 not included in estimate
filterOut.statVis = vis;
filterOut.Phi = Phi;

end