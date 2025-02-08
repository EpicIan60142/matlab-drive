function filterOut = LKF(y, filterParams, filterInit, debug)
% LKF Function that implements a Linearized Kalman Filter for the ASEN 5044
% Statistical Orbit Determination final project
%   - Inputs:
%       - y: Cell array of measurement vectors ordered as follows:
%            {[ [rho^i; rhoDot^i; phi^i; station ID], ...] }
%       - filterParams: Filter parameters structure with the following
%                       fields:
%           - Q: DT Process noise matrix function handle, makes n_wxn_w
%                matrix
%           - R: DT Measurement noise matrix function handle, makes pxp
%                matrix
%       - filterInit: Initial condition structure with the following
%                     fields:
%           - dx0: Initial state perturbation vector for propagation, size 
%                  nx1 array
%           - P0: Initial state covariance matrix, size nxn matrix
%           - xNom: Nominal trajectory for the provided 2 Body Problem,
%                   size nxk matrix
%           - dt: Time step for propagation [sec]
%       - debug: Optional debug vector of booleans defined as follows:
%                [noMeas; noPred]
%           - noMeas: Whether the filter ignores (1) or includes (0)
%                     measurements. Nominally 0 if not specified.
%           - noPred: Whether the filter ignores (1) or includes (0)
%                     measurements. Nominally 0 if not specified.
%   - Outputs:
%       - filterOut: LKF output structure with the following fields:
%           - t_KF: Kalman Filter time vector [sec], size kx1 array
%           - x_KF: Kalman Filter estimated states, size nxk matrix
%           - P_KF: Kalman Filter estimated covariances, size nxnxk matrix
%           - innovation_KF: Kalman Filter innovation vectors, size kx1
%                            cell array of px1 arrays (p can change from
%                            timestep to timestep - multiple stations
%                            possible)
%           - Sk_KF: Kalman Filter innovation covariance matrices, size kx1
%                    cell array of pxp matrices (p can change from timestep
%                    to timestep - multiple stations possible)
%           - sig_KF: 1 sigma standard deviations on state estimates, size
%                     nxk matrix
%
%   Author: Ian Faber, 12/09/2024

if ~iscell(y)
    yCell = {};
    for k = 1:size(y,2)
        meas = [];
        for kk = 1:12
            if ~isnan(y(3*kk-2:3*kk,k))
                meas = [meas, [y(3*kk-2:3*kk,k); kk]];
            else
                meas = meas;
            end
        end
        yCell{k} = meas;
    end
    y = yCell;
end

dT = filterInit.dt;
muEarth = 398600; % km^3/s^2
rEarth = 6378; % km
omegaEarth = (2*pi/86400); % rad/s

    % Propagate stations
numStations = 12; % Number of ground stations tracking the satellite
theta0_s = zeros(1, numStations);
for k = 1:numStations
    theta0_s(k) = (k-1)*(pi/6);
end

    % CT Jacobian function handles
rhoNom = @(x, x_s, y, y_s) sqrt((x-x_s)^2 + (y-y_s)^2);

Atilde = @(x, y) ...
    [
        0                                        1   0                                               0
        muEarth*(2*x^2-y^2)/((x^2+y^2)^(5/2))    0   3*(muEarth*x*y)/((x^2+y^2)^(5/2))               0
        0                                        0   0                                               1
        3*(muEarth*x*y)/((x^2+y^2)^(5/2))        0   muEarth*(2*y^2-x^2)/((x^2+y^2)^(5/2))           0
    ];

Btilde = [
            0   0
            1   0
            0   0
            0   1
         ];

Ctilde = @(x, x_s, xDot, xDot_s, y, y_s, yDot, yDot_s) ... 
    [
         (x-x_s)/rhoNom(x, x_s, y, y_s)                                                        0                                (y-y_s)/rhoNom(x, x_s, y, y_s)                                                        0
         (y-y_s)*((y-y_s)*(xDot-xDot_s) - (x-x_s)*(yDot-yDot_s))/(rhoNom(x, x_s, y, y_s)^3)    (x-x_s)/rhoNom(x, x_s, y, y_s)   (x-x_s)*((x-x_s)*(yDot-yDot_s) - (y-y_s)*(xDot-xDot_s))/(rhoNom(x, x_s, y, y_s)^3)    (y-y_s)/(rhoNom(x, x_s, y, y_s))
         -(y-y_s)/(rhoNom(x, x_s, y, y_s)^2)                                                   0                                (x-x_s)/(rhoNom(x, x_s, y, y_s)^2)                                                    0
    ];

Dtilde = zeros(3,2);

Gamma = [
            0 0
            1 0
            0 0
            0 1
        ];

    % DT Jacobian function handles
F = @(x,y) eye(size(Atilde(x,y))) + dT*Atilde(x,y); % Time varying - depends on x and y!
G = dT*Btilde; % Time invariant
H = Ctilde; % Time varying - depends on x, x_s, xDot, xDot_s, y, y_s, yDot, and yDot_s!
M = Dtilde; % Time invariant
Omega = dT*Gamma; % Time invariant

    % Ground station dynamics function handles
x_stat = @(t,i) rEarth*cos(omegaEarth*t + theta0_s(i));
xDot_stat = @(t,i) -omegaEarth*rEarth*sin(omegaEarth*t + theta0_s(i));
y_stat = @(t,i) rEarth*sin(omegaEarth*t + theta0_s(i));
yDot_stat = @(t,i) omegaEarth*rEarth*cos(omegaEarth*t + theta0_s(i));

    % Helper functions for nominal measurements
rho_stat = @(x, x_s, y, y_s) sqrt((x-x_s)^2 + (y-y_s)^2);
rhoDot_stat = @(x, x_s, xDot, xDot_s, y, y_s, yDot, yDot_s) ...
                ((x-x_s)*(xDot-xDot_s) + (y-y_s)*(yDot-yDot_s))/rho_stat(x, x_s, y, y_s);
phi_stat = @(x, x_s, y, y_s) atan2(y-y_s, x-x_s);

    % Parse debug options
if ~exist("debug", "var")
    noMeas = 0;
    noPred = 0;
else
    noMeas = debug(1);
    noPred = debug(2);
end

    % Pull out filter parameters
Q = filterParams.Q;
R = filterParams.R;

    % Pull out initial conditions, nominal trajectories, and timestep
dx0 = filterInit.dx0;
P0 = filterInit.P0;
xNom = filterInit.xNom;
dt = filterInit.dt;

t_l = 0:dt:14000;

stations_l = [];
for k = 1:length(theta0_s)
    station = struct('id', k, 't', [], 'x', [], 'xDot', [], 'y', [], 'yDot', []);
    stations_l = [stations_l; station];
end

for k = 1:length(theta0_s)
    for kk = 1:length(t_l)
        x_s = x_stat(t_l(kk),k);
        xDot_s = xDot_stat(t_l(kk), k);
        y_s = y_stat(t_l(kk),k);
        yDot_s = yDot_stat(t_l(kk), k);
    
        stations_l(k).t = [stations_l(k).t; t_l(kk)];
        
        stations_l(k).x = [stations_l(k).x; x_s];
        stations_l(k).xDot = [stations_l(k).xDot; xDot_s];
        stations_l(k).y = [stations_l(k).y; y_s];
        stations_l(k).yDot = [stations_l(k).yDot; yDot_s];
    end
end

stations = stations_l;

    % Initialize filter outputs
t_KF = zeros(length(y),1);
x_KF = zeros(size(dx0,1),length(y));
P_KF = zeros(size(P0,1), size(P0,2), length(y));
innovation_KF = cell.empty(length(y),0);
Sk_KF = cell.empty(length(y),0);
sig_KF = zeros(size(x_KF));
IDs_KF = {};
dy_KF = {};
yMeas_KF = {};
yNom_KF = {};

t_KF(1) = 0;
x_KF(:,1) = xNom(:,1) + dx0;
P_KF(:,:,1) = P0;
innovation_KF{1} = zeros(3,1);
Sk_KF{1} = eye(3);
sigComp = [];
for k = 1:size(P_KF, 1)
    sigComp = [sigComp; sqrt(P_KF(k,k,1))];
end
sig_KF(:,1) = sigComp;
IDs_KF{1} = 1;

    % Propagate filter
dx_kPlus = dx0;
P_kPlus = P0;
for k = 1:length(y)-1 % Go to length(y) - 1 because correction step takes data from time k+1
        % Pull out nominal spacecraft state at time k
    x_nom = xNom(1,k);
    xDot_nom = xNom(2,k);
    y_nom = xNom(3,k);
    yDot_nom = xNom(4,k);

    % Time Update / Prediction Step
        % Define dynamics matrices
    F_k = F(x_nom, y_nom);
    Omega_k = Omega;
    Q_k = Q;
    if noPred % Don't include prediction
        dx_kp1Minus = dx_kPlus;
        P_kp1Minus = P_kPlus;
    else
        dx_kp1Minus = F_k*dx_kPlus;
        P_kp1Minus = F_k*P_kPlus*F_k' + Q_k;
    end

    % Measurement Update / Correction Step
    if noMeas % Don't correct if measurements are turned off
        dx_kp1Plus = dx_kp1Minus;
        P_kp1Plus = P_kp1Minus;
    else
            % Define Measurement matrices
        H_kp1 = [];
        R_kp1 = [];
        yMeas = [];
        yNom = [];
        IDcomp = [];
        for kk = 1:size(y{k+1},2) % Account for varying p at each timestep
                % Pull out measurements at time k+1
            yMeas = [yMeas; y{k+1}(1:3,kk)];
            id = y{k+1}(4,kk);

            IDcomp = [IDcomp, id];

                % Pull out nominal spacecraft state at time k+1
            x_nom = xNom(1,k+1);
            xDot_nom = xNom(2,k+1);
            y_nom = xNom(3,k+1);
            yDot_nom = xNom(4,k+1);

                % Pull out station state at time k+1
            x_s = stations(id).x(k+1);
            xDot_s = stations(id).xDot(k+1);
            y_s = stations(id).y(k+1);
            yDot_s = stations(id).yDot(k+1);

                % Create nominal measurement at time k+1
            rho = rho_stat(x_nom, x_s, y_nom, y_s);
            rhoDot = rhoDot_stat(x_nom, x_s, xDot_nom, xDot_s, y_nom, y_s, yDot_nom, yDot_s);
            phi = phi_stat(x_nom, x_s, y_nom, y_s);

            yNom = [yNom; rho; rhoDot; phi];

                % Build up H and R at time k+1
            H_kp1 = [H_kp1; H(x_nom, x_s, xDot_nom, xDot_s, y_nom, y_s, yDot_nom, yDot_s)];
            R_kp1 = blkdiag(R_kp1, R);
        end

            % Make measurement perturbation vector
        dy_kp1 = yMeas - yNom;

        if ~isempty(dy_kp1)
            dy_kp1(3) = wrapToPi(dy_kp1(3));
            if size(dy_kp1,1) > 3
                dy_kp1(6) = wrapToPi(dy_kp1(6));
            end
        end

        % for points = [55, 102, 149, 196, 243, 290, 346, 391, 437, 485, 533]
        %     if k == points
        %         fprintf("Timestep: %.0f, yMeas: [%.04f; %.04f; %.04f], yNom: [%.04f; %.04f; %.04f], dy: [%.04f; %.04f; %.04f]\n", k, yMeas, yNom, dy_kp1)
        %     end
        % end

            % Update state estimate
        if isempty(y{k+1})
            R_kp1 = R;
            dx_kp1Plus = dx_kp1Minus;
            P_kp1Plus = P_kp1Minus;
        else
                % Define Kalman Gain matrix
            K_kp1 = P_kp1Minus*H_kp1'*(H_kp1*P_kp1Minus*H_kp1' + R_kp1)^-1;

                % Get next state and covariance estimates
            dx_kp1Plus = dx_kp1Minus + K_kp1*(dy_kp1 - H_kp1*dx_kp1Minus);
            P_kp1Plus = (eye(size(K_kp1,1)) - K_kp1*H_kp1)*P_kp1Minus;
        end
    end

    % Save outputs
    t_KF(k+1) = k*dt;
    x_KF(:,k+1) = xNom(:,k+1) + dx_kp1Plus;
    P_KF(:,:,k+1) = P_kp1Plus;
    if isempty(y{k+1}) || noMeas
        innovation_KF{k+1} = [];%zeros(3,1);
        Sk_KF{k+1} = R_kp1;
    else
        innovation_KF{k+1} = dy_kp1 - H_kp1*dx_kp1Minus;
        Sk_KF{k+1} = H_kp1*P_kp1Minus*H_kp1' + R_kp1;
    end
    
    sigComp = [];
    for kk = 1:size(P_KF, 1)
        sigComp = [sigComp; sqrt(P_KF(kk,kk,k+1))];
    end
    sig_KF(:,k+1) = sigComp;
    IDs_KF{k+1} = IDcomp;
    dy_KF{k+1} = dy_kp1;
    yMeas_KF{k+1} = yMeas;
    yNom_KF{k+1} = yNom;

    % Update for next run
    dx_kPlus = dx_kp1Plus;
    P_kPlus = P_kp1Plus;
end

    % Assign function outputs
filterOut.t_KF = t_KF;
filterOut.x_KF = x_KF;
filterOut.P_KF = P_KF;
filterOut.innovation_KF = innovation_KF;
filterOut.Sk_KF = Sk_KF;
filterOut.sig_KF = sig_KF;
filterOut.IDs_KF = IDs_KF;
filterOut.stations = stations;
filterOut.dy_KF = dy_KF;
filterOut.yMeas_KF = yMeas_KF;
filterOut.yNom_KF = yNom_KF;

end