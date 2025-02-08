function filterOut = KalmanFilter(y, F, Q, H, R, x0, P0, dt, debug)
% Function that implements a Kalman Filter on some set of existing data
%   Inputs:
%       - y: Vector of measurements
%       - F: DT Dynamics matrix, size nxn
%       - Q: Kalman Filter DT process noise intensity matrix, size nxn
%       - H: DT Output matrix, size pxn
%       - R: Kalman Filter DT measurement noise intensity matrix, size pxp
%       - x0: Initial state for propagation, size nx1
%       - P0: Initial state covariance, size nxn
%       - dt: Time step for propagation [sec]
%       - debug: Debug vector of booleans as follows:
%           [noMeas; noPred]
%           - noMeas: Whether the filter ignores (1) or includes (0) 
%                     measurements
%           - noPred: Whether the filter ignores (1) or includes (0)
%                     predictions
%
%   Outputs:
%       - filterOut: Kalman Filter output structure with the following
%                    fields:
%           - t_KF: Kalman Filter time vector [sec]
%           - x_KF: Kalman Filter estimated states
%           - P_KF: Kalman Filter estimated covariances
%           - innovation_KF: Kalman Filter innovation vectors
%           - Sk_KF: Kalman Filter innovation covariance matrices
%           - sig_KF: 1 sigma standard deviations on state estimates
%
%   By: Ian Faber, 11/18/2024
%

noMeas = debug(1);
noPred = debug(2);

t_KF = [];
x_KF = [];
P_KF = zeros(size(F,1), size(F,2), size(y,2));
innovation_KF = [];
Sk_KF = zeros(size(H,1), size(H,1), size(y,2));
sig_KF = [];

x_kPlus = x0;
P_kPlus = P0;
for k = 1:size(y,2)
    % Matrix assignment
    H_kp1 = H;
    R_kp1 = R;

    % Time Update / Prediction Step
    if noPred % Don't include prediction
        x_kp1Minus = x_kPlus;
        P_kp1Minus = P_kPlus;
    else
        x_kp1Minus = F*x_kPlus;
        P_kp1Minus = F*P_kPlus*F' + Q;
    end
    K_kp1 = P_kp1Minus*H_kp1'*(H_kp1*P_kp1Minus*H_kp1' + R_kp1)^-1;

    % Measurement Update / Correction Step
    if noMeas % Don't include measurements
        x_kp1Plus = x_kp1Minus;
        P_kp1Plus = P_kp1Minus;
    else
        x_kp1Plus = x_kp1Minus + K_kp1*(y(:,k) - H_kp1*x_kp1Minus); % y starts at y_1, not y_0
        P_kp1Plus = (eye(size(F))-K_kp1*H_kp1)*P_kp1Minus;
    end

    % Save outputs
    t_KF = [t_KF; k*dt];
    x_KF = [x_KF, x_kp1Plus];
    P_KF(:,:,k) = P_kp1Plus;
    innovation_KF = [innovation_KF, y(:,k) - H_kp1*x_kp1Minus];
    Sk_KF(:,:,k) = H_kp1*P_kp1Minus*H' + R_kp1;
    sigComp = [];
    for kk = 1:size(P_KF, 1)
        sigComp = [sigComp; sqrt(P_KF(kk,kk,k))];
    end
    sig_KF = [sig_KF, sigComp];

    % Update for next run
    x_kPlus = x_kp1Plus;
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