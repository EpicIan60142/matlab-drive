function Qw = makeQ_DMC(B,Qu,t,t0)
% Function that calculates the 9x9 Q matrix for models including dynamic
% model compensation
%   - Inputs:
%       - B: DMC time constant matrix organized as a 3x3 matrix:
%            B = diag([tau_x^-1, tau_y^-1, tau_z^-1])
%       - Qu: Process noise covariance matrix organized as a 3x3 matrix:
%             Qu = diag([sigma_x^2, sigma_y^2, sigma_z^2])
%       - t: Current time to generate the Q matrix at
%       - t0: Initial time to reference for generating Q - this is
%             generally t_im1 in the LKF/EKF
%
%   By: Ian Faber, 02/24/2025
%

Qw = [];
for k = 1:size(B,1)
    beta_i = B(k,k);
    sigSqrd_i = Qu(k,k);

        % Qw(r,r)
    Q11 = sigSqrd_i*((1/(3*beta_i^2))*(t-t0)^3 - (1/beta_i^3)*(t-t0)^2 + (1/beta_i^4)*(t-t0) - (2/beta_i^4)*exp(-beta_i*(t-t0))*(t-t0) + (1/(2*beta_i^5))*(1-exp(-2*beta_i*(t-t0))));
        % Qw(r,v) = Qw(v,r)
    Q21 = sigSqrd_i*((1/(2*beta_i^2))*(t-t0)^2 - (1/beta_i^3)*(t-t0) + (1/beta_i^3)*exp(-beta_i*(t-t0))*(t-t0) + (1/beta_i^4)*(1-exp(-beta_i*(t-t0))) - (1/(2*beta_i^4))*(1-exp(-2*beta_i*(t-t0))));
    Q12 = Q21;
        % Q2(r,w) = Qw(w,r)
    Q31 = sigSqrd_i*((1/(2*beta_i^3))*(1-exp(-2*beta_i*(t-t0))) - (1/beta_i^2)*exp(-beta_i*(t-t0))*(t-t0));
    Q13 = Q31;
        % Qw(v,v)
    Q22 = sigSqrd_i*((1/beta_i^2)*(t-t0) - (2/beta_i^3)*(1-exp(-beta_i*(t-t0))) + (1/(2*beta_i^3))*(1-exp(-2*beta_i*(t-t0))));
        % Qw(v,w) = Qw(w,v)
    Q23 = sigSqrd_i*((1/(2*beta_i^2))*(1+exp(-2*beta_i*(t-t0))) - (1/beta_i^2)*exp(-beta_i*(t-t0)));
    Q32 = Q23;
        % Qw(w,w)
    Q33 = (sigSqrd_i/(2*beta_i))*(1-exp(-2*beta_i*(t-t0)));

    Q_i = [
            Q11, Q12, Q13;
            Q21, Q22, Q23;
            Q31, Q32, Q33
          ];
    
    Qw = blkdiag(Qw, Q_i);
end

end

