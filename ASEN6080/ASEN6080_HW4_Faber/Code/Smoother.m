function smoothOut = Smoother(filterOut)
% Function that implements a sequential smoothing algorithm for Stat OD
% problems
%   - Inputs:
%       - filterOut: Output structure from an LKF run, as defined in
%                    LKF_SNC.m
%   - Outputs:
%       - smoothOut: Smoother output structure with the following fields:
%           - tSmoothed: Smoothed estimate time
%           - xSmoothed: Smoothed deviation estimates:
%                        [xSmoothed_1, xSmoothed_2, ..., xSmoothed_k]
%           - PSmoothed: Smoothed state covariance estimates:
%                        [{PSmoothed_1}, {PSmoothed_2}, ..., {PSmoothed_k}]
%           - XSmoothed: Smoothed full state estimate, calculated as 
%                        XSmoothed = Xstar + xSmoothed:
%                        [XSmoothed_1, XSmoothed_2, ..., XSmoothed_k]
%
%   By: Ian Faber, 03/03/2025
%
    % Process filterOut structure
t = filterOut.t;
xEst = filterOut.xEst;
PBarEst = filterOut.PBarEst;
PEst = filterOut.PEst;
XStar = filterOut.XStar;
Phi = filterOut.Phi;

    % Initialize outputs
tSmoothed = [];
xSmoothed = [];
PSmoothed = [];
XSmoothed = [];

    % Find value of l
l = size(xEst, 2);

    % Run smoother
x_kp1l = xEst(:,l); % Deviation at t_kp1 (t_l) given measurements at t_l
P_kp1l = PEst{l}; % Covariance at t_kp1 (t_l) given measurements at t_l
for k = (l-1):-1:1
        % Pull out filter estimates
    x_kk = xEst(:,k); % Deviation at t_k given measurements at t_k
    P_kp1k = PBarEst{k+1}; % Covariance at t_kp1 given measurements at t_k
    P_kk = PEst{k}; % Covariance at t_k given measurements at t_k
    tSmoothed = [tSmoothed; t(k)];

        % Pull out Phi(t_kp1, t_k)
    Phi_kp1 = Phi{k+1};

        % Smooth estimate
    % S_k = P_kk*Phi_kp1'*(Phi_kp1*P_kk*Phi_kp1' + Q_kp1)^-1;
    S_k = P_kk*Phi_kp1'*(P_kp1k)^-1;
    x_kl = x_kk + S_k*(x_kp1l - Phi_kp1*x_kk);
        
        % Smooth covariance
    P_kl = P_kk + S_k*(P_kp1l - P_kp1k)*S_k';
        
        % Save outputs
    xSmoothed = [xSmoothed, x_kl];
    PSmoothed = [PSmoothed; {P_kl}];
    XSmoothed = [XSmoothed, XStar(:,k) + x_kl];

        % Update for next run
    x_kp1l = x_kl;
    P_kp1l = P_kl;
    
    
end

    % Assign smoother outputs - need to flip them to be time ascending!
smoothOut.tSmoothed = flipud(tSmoothed);
smoothOut.xSmoothed = fliplr(xSmoothed);
smoothOut.PSmoothed = flipud(PSmoothed);
smoothOut.XSmoothed = fliplr(XSmoothed);

end