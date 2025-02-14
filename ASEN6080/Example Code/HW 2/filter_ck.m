function [state_out, Covariance, residuals, resid_pf, P_pf, X_pf] = filter_ck(x_init, Pbar_init, obs, obsW, options )
%
% ---------------------------------------------------------------------
% 
% Description:
%
%   Implementation of a generic Conventional Kalman Filter.  Intended to be
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
%       options.joseph        - flag to indicate Joseph Formulation
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
%       Pbar_init - initial covariance matrix
%       obs       - observations structure.  requires array of times and
%                       observations
%       obsW      - observations weighting matrix
%       options   - options structure as previously described
% 
% Outputs:
% 
%       state_out  - CKF determined state at the final time
%       Covariance - final post-fit covariance
%       residuals  - pre-fit residuals
%       resid_pf   - post-fit residuals
%       P_pf       - reshaped covariance at each observation time
%       X_pf       - Post-fit state solution
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
%   01feb06     Brandon A. Jones      original version, adapted from old
%                                       code
%   22feb06     Brandon A. Jones      added ability to pass arguments to
%                                       miscellaneous functions.
%   2025        Jay McMahon           various modifications over the years

num_states   = length(x_init);
Phi0         = eye(num_states);
Phi0_reshape = reshape(Phi0,num_states*num_states,1);
x0           = x_init;
Pbar0        = Pbar_init;
xbar0        = zeros(num_states,1);
R            = inv(obsW);

abs_tol     = ones(1,num_states+num_states*num_states).*options.tol;
ode_options = odeset('RelTol',options.tol,'AbsTol',abs_tol);
CovarianceI = eye(num_states);

status_chars = ['-' '\' '|' '/'];
fprintf('%c',status_chars(1));

for iterations = 1:options.num_iters
    for j = 1:length(obs.times(:))
        
        fprintf('\b%c',status_chars(mod(j,4)+1));

        if j==1

            x_hat      = xbar0;
            Covariance = Pbar0;
            prev_time  = 0;
            PhiTotal   = Phi0;
            residuals  = zeros(length(obs.times(:)),length(obs.obs(1,:)));
            resid_pf   = zeros(length(obs.times(:)),length(obs.obs(1,:)));
            P_pf       = zeros(length(obs.times(:)),num_states^2);
            X_pf       = zeros(length(obs.times(:)),num_states);

            X_Int(1,:) = [x0; Phi0_reshape ];
            Time_Int   = 0;
        else
            [Time_Int,X_Int] = ode45(options.integ_fcn,[prev_time obs.times(j)],...
                [X_Ref; Phi0_reshape], ode_options, options.integ_args{:});
            prev_time = obs.times(j);
        end

        prev_time = obs.times(j);
        Phi = reshape(X_Int(end,(num_states+1):end),num_states,num_states);

        PhiTotal = Phi*PhiTotal;

        %  Time Update
        xbar = Phi*x_hat;
        Pbar = Phi*Covariance*Phi';

        C = feval(options.obs_fcn, obs, X_Int(end,:), j, Time_Int(end), options.extra_args );
        OminusC = (obs.obs(j,:)-C)';
        residuals(j,:) = OminusC;
        

        Htilde = feval(options.H_fcn, X_Int(end,:), obs, obs.times(j), j, options.extra_args );
        KalScratch = Pbar*Htilde';

        KalGain = KalScratch/(Htilde*KalScratch + R);
        
        %  Meas Update
        x_hat = xbar + KalGain*(OminusC - Htilde*xbar);

        if options.joseph == 1
            Covariance = (CovarianceI-KalGain*Htilde)*Pbar*(CovarianceI-KalGain*Htilde)' ...
                + KalGain*R*KalGain';
        else
            Covariance = (CovarianceI - KalGain*Htilde)*Pbar;
        end
        
        X_Ref = X_Int(end,1:num_states)';
        
        P_pf(j,:) = reshape(Covariance,1,num_states^2);
        
        % Post-fit residuals
        resid_pf(j,:) = OminusC - Htilde*x_hat;
        
        % Post-fit state solution
        X_pf(j,:) = (X_Ref + x_hat)';

    end

    x_hat0 = PhiTotal\x_hat;
    
    x0    = x0 + x_hat0;
    xbar0 = xbar0 - x_hat0;

    clear X_Int

end

residuals = [obs.times residuals];

state_out = X_pf(end,:)';

fprintf('\b');
