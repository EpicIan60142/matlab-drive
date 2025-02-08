function filterOut = LKF(y, filterParams, filterInit, debug)
% LKF Function that implements a Linearized Kalman Filter for the ASEN 5044
% Statistical Orbit Determination final project
%   - Inputs:
%       - y: Cell array of measurement vectors ordered as follows:
%            {[ [rho^i; rhoDot^i; phi^i; station ID], ...] }
%       - filterParams: LKF input structure with the following fields:
%           - F: DT Dynamics Matrix function handle, creates size nxn
%                matrix
%           - Q: DT Process Noise Intensity Matrix function handle, creates
%                size n_wxn_w matrix
%           - Omega: DT Process Noise Dyanmics Matrix function handle,
%                    creates size nxn_w matrix
%           - H: DT Output Matrix function handle, creates size pxn matrix
%           - R: DT Measurement Noise Intensity Matrix function handle,
%                creates size pxp matrix
%           - dt: Time step for propagation [sec], scalar
%       - filterInit: Initial condition structure with the following
%                     fields:
%           - dx0: Initial state perturbation vector for propagation, size 
%                  nx1 array
%           - P0: Initial state covariance matrix, size nxn matrix
%           - xNom: Nominal trajectory for the provided 2 Body Problem,
%                   size nxk matrix
%           - stations: 12x1 structure of station states at each time step 
%                       with the following fields for each station:
%               - id: Station ID, scalar from 1-12
%               - t: time vector [sec], size kx1 array
%               - x: X coordinate of station i, size kx1 array
%               - xDot: X velocity of station i, size kx1 array
%               - y: Y coordinate of station i, size kx1 array
%               - yDot: Y velocity of station i, size kx1 array
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
F = filterParams.F;
Q = filterParams.Q;
Omega = filterParams.Omega;
H = filterParams.H;
R = filterParams.R;
dt = filterParams.dt;

    % Pull out initial conditions, nominal trajectories, and stations
dx0 = filterInit.dx0;
P0 = filterInit.P0;
xNom = filterInit.xNom;
stations = filterInit.stations;

    % Initialize filter outputs
t_KF = zeros(length(y),1);
x_KF = zeros(size(dx0,1),length(y));
P_KF = zeros(size(P0,1), size(P0,2), length(y));
innovation_KF = cell.empty(length(y),0);
Sk_KF = cell.empty(length(y),0);
sig_KF = zeros(size(x_KF));

t_KF(1) = 0;
x_KF(:,1) = xNom(:,1) + dx0;
P_KF(:,:,1) = P0;
sigComp = [];
for k = 1:size(P_KF, 1)
    sigComp = [sigComp; sqrt(P_KF(k,k,1))];
end
sig_KF(:,1) = sigComp;

    % Propagate filter
dx_kPlus = dx0;
P_kPlus = P0;
for k = 1:length(y)-1 % Go to length(y) - 1 because correction step takes data from time k+1
        % Pull out nominal spacecraft state at time k
    x_sat = xNom(1,k);
    xDot_sat = xNom(2,k);
    y_sat = xNom(3,k);
    yDot_sat = xNom(4,k);

    % Time Update / Prediction Step
        % Define dynamics matrices
    F_k = F(x_sat, y_sat);
    Omega_k = Omega;
    Q_k = Q;
    if noPred % Don't include prediction
        dx_kp1Minus = dx_kPlus;
        P_kp1Minus = P_kPlus;
    else
        dx_kp1Minus = F_k*dx_kPlus;
        P_kp1Minus = F_k*P_kPlus*F_k' + Omega_k*Q_k*Omega_k';
    end

    % Measurement Update / Correction Step
    if noMeas % Don't correct if measurements are turned off
        dx_kp1Plus = dx_kp1Minus;
        P_kp1Plus = P_kp1Minus;
    elseif isempty(y{k+1}) % Don't correct if measurement is empty
        dx_kp1Plus = dx_kp1Minus;
        P_kp1Plus = P_kp1Minus;
    else
            % Define Measurement matrices
        H_kp1 = [];
        R_kp1 = [];
        yMeas = [];
        yNom = [];
        for kk = 1:size(y{k+1},2) % Account for varying p at each timestep
                % Pull out measurements at time k+1
            yMeas = [yMeas; y{k+1}(1:3,kk)];
            id = y{k+1}(4,kk);

                % Pull out nominal spacecraft state at time k+1
            x_sat = xNom(1,k+1);
            xDot_sat = xNom(2,k+1);
            y_sat = xNom(3,k+1);
            yDot_sat = xNom(4,k+1);

                % Pull out station state at time k+1
            x_s = stations(id).x(k+1);
            xDot_s = stations(id).xDot(k+1);
            y_s = stations(id).y(k+1);
            yDot_s = stations(id).yDot(k);

                % Create nominal measurement at time k+1
            yNom = [yNom; rho_stat(x_sat, x_s, y_sat, y_s); rhoDot_stat(x_sat, x_s, xDot_sat, xDot_s, y_sat, y_s, yDot_sat, yDot_s); phi_stat(x_sat, x_s, y_sat, y_s)];

                % Build up H and R at time k+1
            H_kp1 = [H_kp1; H(x_sat, x_s, xDot_sat, xDot_s, y_sat, y_s, yDot_sat, yDot_s)];
            R_kp1 = blkdiag(R_kp1, R);
        end
            % Make measurement perturbation vector
        dy_kp1 = yMeas - yNom;

            % Define Kalman Gain matrix
        K_kp1 = P_kp1Minus*H_kp1'*(H_kp1*P_kp1Minus*H_kp1' + R_kp1)^-1;

            % Update state estimate
        dx_kp1Plus = dx_kp1Minus + K_kp1*(dy_kp1 - H_kp1*dx_kp1Minus);
        P_kp1Plus = (eye(size(K_kp1,1)) - K_kp1*H_kp1)*P_kp1Minus;
    end

    % Save outputs
    t_KF(k+1) = k*dt;
    x_KF(:,k+1) = xNom(:,k+1) + dx_kp1Plus;
    P_KF(:,:,k+1) = P_kp1Plus;
    if isempty(y{k+1}) || noMeas
        innovation_KF{k+1} = [];
        Sk_KF{k+1} = [];
    else
        innovation_KF{k+1} = yMeas - H_kp1*dx_kp1Minus;
        Sk_KF{k+1} = H_kp1*P_kp1Minus*H_kp1' + R_kp1;
    end
    
    sigComp = [];
    for kk = 1:size(P_KF, 1)
        sigComp = [sigComp; sqrt(P_KF(kk,kk,k+1))];
    end
    sig_KF(:,k+1) = sigComp;

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

end