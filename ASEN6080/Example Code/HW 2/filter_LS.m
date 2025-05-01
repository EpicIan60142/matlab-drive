function [state_out, Covar, resids, resid_pf, xhat] = filter_LS(x_init, xbar, Pbar_init, obs, obsW, options)
%
% ---------------------------------------------------------------------
% 
% Description:
%
%   Implementation of a generic Least Squares Filter.  Intended to be
%   a plug-and-play interface for quick reconfiguration of the code.
%
%   Implementation is kept generic through the use of an options structure.
%   Currently, there are several variables that must be defined in the
%   structure for proper execution of the filter.  They are:
%
%       options.tol           - integration tolerance for ODE45
%       options.num_iters     - number of iterations of the CKF
%       options.integ_fcn     - function that handles state integration 
%       options.H_fcn         - function that generates the H matrix
%       options.obs_fcn       - function that computes observations
%       options.extra_args    - structure containing extra arguments for
%                               various functions.  this structure must 
%                               include at least one variable, even if it 
%                               is a dummy.
%       options.integ_args    - cell array of extra variables to pass to
%                               the integration function.  must exist, 
%                               but can be empty.
%
%   Any options specific to the function handles provided in the options
%   structure can be included as well.
% 
% Inputs:
%
%       x_init    - initial state vector
%       xbar      - initial state deviation
%       Pbar_init - initial covariance matrix
%       obs       - observations structure.  requires array of times and
%                       observations
%       obsW      - observations weighting matrix
%       options   - options structure as previously described
% 
% Outputs:
% 
%       state_out - LS determined state at the epoch initial time
%       Covar     - variance-covariance matrix of the estimate
%       resids    - pre-fit measurement residuals
%       resid_pf  - post-fit measurement residuals
%       xhat      - final state deviation estimate at epoch
%
% Assumptions/References:
%
%   Based on algorithm in "Statistical Orbit Determination" by Tapley,
%   Schutz and Born, p. 203.
%
%   Assumes the options structure has been correctly populated as described
%   previously.
%
% Dependencies:
%
%   ode45.m
%   any functions provided in the options structure
%
% Modification History:
% 
%   01feb06     Brandon A. Jones      original version
%   2025        Jay McMahon           various modifications over the years
%        

%  Organize the inputs for use in the filter.
num_states    = length(x_init);
Phi0          = eye(num_states);
x0            = [x_init; reshape(Phi0,num_states*num_states,1)];
Pbar0         = Pbar_init;

%  Set the time and filter parameters required for execution
time          = options.start_time:options.inc_time:options.end_time;
abs_tol       = ones(1,length(x0)).*options.tol;
converge_crit = options.conv_crit;
ode_options   = odeset( 'RelTol', options.tol, 'AbsTol', abs_tol);

j = 0;
xhat_mag = 10*converge_crit;

while xhat_mag > converge_crit
    
    j = j+1;

    [t,x]       = ode45( options.integ_fcn, time, x0, ode_options, options.integ_args{:});

    
    PbarR = chol(Pbar0);
    PbarR_inv = inv(PbarR);
    
    Pbar_inv = PbarR_inv*PbarR_inv';
    
    L = Pbar_inv;
    N = Pbar_inv*xbar;

    OminusC = zeros(length(obs.times),2);
    for i = 1:length(obs.times)
        index = find(obs.times(i)== t );
        if length(index) < 1
            continue
        end
        Phi = reshape(x(index,(num_states+1):end),num_states,num_states);
        C = feval(options.obs_fcn, obs, x(index,1:num_states), i, t(index), options.extra_args );
        OminusC(i,:) = obs.obs(i,:)-C;
        Htilde = feval(options.H_fcn, x(index,1:num_states), obs, t(index), i, options.extra_args );
        H = Htilde*Phi;
        HtW = H'*obsW;
        L = L + HtW*H;
        N = N + HtW*OminusC(i,:)';
    end
    
    LR = chol(L);
    
    LR_inv = inv(LR);
    Covar = LR_inv*LR_inv';

    xhat = Covar*N;

    xhat_mag = norm(xhat);

    x0(1:num_states) = x0(1:num_states) + xhat;

    % Pre-fits RMS
    RMS_Rho = sqrt(sum(OminusC(:,1).*OminusC(:,1))/length(OminusC(:,1)));
    RMS_RhoDot = sqrt(sum(OminusC(:,2).*OminusC(:,2))/length(OminusC(:,2)));

    % Compute linearized post-fits
    resid_pf = zeros(length(obs.times),2);
    for i = 1:length(obs.times)
        index = find(obs.times(i)== t );
        if length(index) < 1
            continue
        end
        Phi = reshape(x(index,(num_states+1):end),num_states,num_states);
        Htilde = feval(options.H_fcn, x(index,1:num_states), obs, t(index), i, options.extra_args );
        H = Htilde*Phi;
        resid_pf(i,:) = OminusC(i,:) - (H*xhat)';
    end
    RMS_Rho_pf = sqrt(sum(resid_pf(:,1).*resid_pf(:,1))/length(resid_pf(:,1)));
    RMS_RhoDot_pf = sqrt(sum(resid_pf(:,2).*resid_pf(:,2))/length(resid_pf(:,2)));
    
    fprintf('--- Iteration %d info: ---\n', j );
    fprintf('Pre-fit residual RMS values\n');
    fprintf('Rho    = %g meters\n', RMS_Rho );
    fprintf('RhoDot = %g meters\n', RMS_RhoDot );
    fprintf('Post-fit residual RMS values\n');
    fprintf('Rho    = %g meters\n', RMS_Rho_pf );
    fprintf('RhoDot = %g meters\n', RMS_RhoDot_pf );
    fprintf('x_hat magnitude:  %g\n', xhat_mag );
    fprintf('\n---------------------------------------------\n\n');

    % Make plots per iteration
    resids = [obs.times OminusC];
    sta1Ind = find(obs.station==1);
    sta2Ind = find(obs.station==2);
    sta3Ind = find(obs.station==3);

    figure();
    set(gcf,'Color','w');
    subplot(2,1,1);
    plot(resids(sta1Ind,1)./(3600),resids(sta1Ind,2),'.',resids(sta2Ind,1)./(3600),resids(sta2Ind,2),'.',resids(sta3Ind,1)./(3600),resids(sta3Ind,2),'.'); % add different colors for different stations
    title(['Pre-fit Residuals: Iteration ' num2str(j)])
    ylabel('Range [km]');
    subplot(2,1,2);
    plot(resids(sta1Ind,1)./(3600),resids(sta1Ind,3),'.',resids(sta2Ind,1)./(3600),resids(sta2Ind,3),'.',resids(sta3Ind,1)./(3600),resids(sta3Ind,3),'.');
    ylabel('Range-Rate [km/s]')
    xlabel('Time (hours)')

    figure();
    set(gcf,'Color','w');
    subplot(2,1,1);
    plot(resids(sta1Ind,1)./(3600),resid_pf(sta1Ind,1),'.',resids(sta2Ind,1)./(3600),resid_pf(sta2Ind,1),'.',resids(sta3Ind,1)./(3600),resid_pf(sta3Ind,1),'.'); % add different colors for different stations
    title(['Post-fit Residuals: Iteration ' num2str(j)])
    ylabel('Range [km]');
    subplot(2,1,2);
    plot(resids(sta1Ind,1)./(3600),resid_pf(sta1Ind,2),'.',resids(sta2Ind,1)./(3600),resid_pf(sta2Ind,2),'.',resids(sta3Ind,1)./(3600),resid_pf(sta3Ind,2),'.');
    ylabel('Range-Rate [km/s]')
    xlabel('Time (hours)')
    
    % Update prior deviation estimate
    xbar = xbar - xhat;
    
    if j >= options.max_iterations
        break;
    end

end

state_out = x0(1:num_states);

resids = [obs.times OminusC];

