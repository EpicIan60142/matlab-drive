tic

%% Quick User Inputs (basic case changing)
observationFile = 'HW2_obs.txt';
useFilter = 3; % 1 for batch, 2 for CKF, 3 for EKF
dx = zeros(6,1);
sigma_r = 1;    % km
sigma_v = 1e-3; % km/s
dx = [1; 1; 1; 1e-3; 1e-3; 1e-3];


%% Set up problem
% a priori state, covariance
a       = 10000.0; % km
e       = 0.001;
MU         = 3.986004415E5; % km^3/s^2
[r, v]  = kep2cart(a,e,40.0*pi/180,80.0*pi/180,40.0*pi/180,0,MU);
num_states = 6;
x_init = [r'; v'] + dx;
xbar0 = zeros(num_states,1); 
Pbar_init  = [sigma_r^2.*eye(3) zeros(3); zeros(3) sigma_v^2.*eye(3)];

% Dynamics Models
options.integ_fcn     = @two_body_wPhi;
J2         = 0.0010826269;
REarth     = 6378.0; % km
tdot       = 7.2921158553E-5; % rad/s
theta0     = 122*pi/180.0; % rad

options.integ_args = {REarth,MU,J2};

% Measurement Models
options.obs_fcn       = @epoch_range_rangerate;
options.H_fcn         = @Htilde_5070;

%  Set the observation variance-covariance matrix
sigma_rho    = 1e-3; % km
sigma_rhodot = 1e-6; % km/s
obsW         = [1/(sigma_rho*sigma_rho) 0; 0 1/(sigma_rhodot*sigma_rhodot)];

% Station locations
[xs,ys,zs] = sph2cart(148.981944*pi/180, -35.398333*pi/180, 6378);
Rs1 = [xs;ys;zs];
[xs,ys,zs] = sph2cart(355.749444*pi/180, 40.427222*pi/180, 6378);
Rs2 = [xs;ys;zs];
[xs,ys,zs] = sph2cart(243.205*pi/180, 35.247164*pi/180, 6378);
Rs3 = [xs;ys;zs];

%  Set the extra arguments required for evaluating the Htilde and computed
%  observations
options.extra_args.tdot       = tdot;
options.extra_args.theta0     = theta0;
options.extra_args.S1         = Rs1;
options.extra_args.S2         = Rs2;
options.extra_args.S3         = Rs3;

% Load data
ObsData = importdata(observationFile);
obs.times   = ObsData(:,1);
obs.station = ObsData(:,2);
obs.obs     = ObsData(:,3:4);


%% Set up filter

% Common Options
options.tol           = 1e-10;

% Batch Options
options.start_time     = 0;
options.end_time       = max(obs.times);
options.inc_time       = 10;
options.conv_crit      = 0.001; %0.3;
options.max_iterations = 10;

% CKF Options
options.num_iters     = 3; % do this many iterations of KF
options.joseph        = 1;

% EKF Options
options.num_obs_switch = 100;

%% Run the filter
switch useFilter
    
    case 1 % batch
        [state_out, covar, resids, resids_pf, xhat] = filter_LS(x_init, xbar0, Pbar_init, obs, obsW, options);
        
    case 2 % CKF
        [state_out, covar, resids, resids_pf, P_pf, X_pf] = filter_ck(x_init, Pbar_init, obs, obsW, options );
        
    case 3 % EKF
        [state_out, covar, resids, resids_pf, P_pf, X_pf] = filter_ekf(x_init, Pbar_init, obs, obsW, options );
        
end



%% Outputs and analysis
% Choose what to run
%   In some cases, won't need to propagate or load truth trajectories


% Load truth data
haveTruth = 1;
truthTrajFile = 'HW2_truth_J2.txt';
if haveTruth == 1
    Xtrue = importdata(truthTrajFile);
end

% Construct best fit trajectory from outputs
abs_tol       = ones(1,length(state_out)+36).*options.tol;
ode_options   = odeset( 'RelTol', options.tol, 'AbsTol', abs_tol);
switch useFilter
    
    case 1 % batch
        
        % solution at epoch, propagate forward
        [~,X_bestfit]       = ode45( options.integ_fcn, [options.start_time:options.inc_time:options.end_time], [state_out; reshape(eye(6),36,1)], ode_options, options.integ_args{:});
        % Propagate covariance
        Phist = zeros(length(X_bestfit),num_states^2);
        for ii = 1:length(X_bestfit)
            Phi = reshape(X_bestfit(ii,(num_states+1):end),num_states,num_states);
            Pmapped = Phi*covar*Phi';
            Phist(ii,:) = reshape(Pmapped,1,num_states^2);
        end
        
    otherwise % CKF or EKF without smoother
        
        % solution at final time; propagate backwards
        [~,Xprop]       = ode45( options.integ_fcn, fliplr([options.start_time:options.inc_time:options.end_time]), [state_out; reshape(eye(6),36,1)], ode_options, options.integ_args{:});
        X_bestfit = flipud(Xprop);
        % Propagate covariance
        Phist = zeros(length(X_bestfit),num_states^2);
        for ii = 1:length(Xprop)
            Phi = reshape(Xprop(ii,(num_states+1):end),num_states,num_states);
            Pmapped = Phi*covar*Phi';
            Phist(ii,:) = reshape(Pmapped,1,num_states^2);
        end
        Phist = flipud(Phist);
            
end

% Run analysis and plots
% Can adjust starting index to only compute steady state RMS
RMS_zero_ind = 1;

RMS_rho_pf = sqrt(sum(resids_pf(RMS_zero_ind:end,1).*resids_pf(RMS_zero_ind:end,1))/length(resids_pf(RMS_zero_ind:end,1)));
RMS_rhodot_pf = sqrt(sum(resids_pf(RMS_zero_ind:end,2).*resids_pf(RMS_zero_ind:end,2))/length(resids_pf(RMS_zero_ind:end,2)));

fprintf('\nPost-fit residuals\n');
fprintf('Rho    = %g km\n', RMS_rho_pf );
fprintf('RhoDot = %g km/s\n', RMS_rhodot_pf );

fprintf('--------------------------------------------------------\n');

%  Generate profiling plots

sta1Ind = find(obs.station==1);
sta2Ind = find(obs.station==2);
sta3Ind = find(obs.station==3);

figure();
set(gcf,'Color','w');
subplot(2,1,1);
plot(resids(sta1Ind,1)./(3600),resids(sta1Ind,2),'.',resids(sta2Ind,1)./(3600),resids(sta2Ind,2),'.',resids(sta3Ind,1)./(3600),resids(sta3Ind,2),'.'); % add different colors for different stations
title('Pre-fit Residuals')
ylabel('Range [km]');
subplot(2,1,2);
plot(resids(sta1Ind,1)./(3600),resids(sta1Ind,3),'.',resids(sta2Ind,1)./(3600),resids(sta2Ind,3),'.',resids(sta3Ind,1)./(3600),resids(sta3Ind,3),'.');
ylabel('Range-Rate [km/s]')
xlabel('Time (hours)')

figure();
set(gcf,'Color','w');
subplot(2,1,1);
plot(resids(sta1Ind,1)./(3600),resids_pf(sta1Ind,1),'.',resids(sta2Ind,1)./(3600),resids_pf(sta2Ind,1),'.',resids(sta3Ind,1)./(3600),resids_pf(sta3Ind,1),'.'); % add different colors for different stations
title('Post-fit Residuals')
ylabel('Range [km]');
subplot(2,1,2);
plot(resids(sta1Ind,1)./(3600),resids_pf(sta1Ind,2),'.',resids(sta2Ind,1)./(3600),resids_pf(sta2Ind,2),'.',resids(sta3Ind,1)./(3600),resids_pf(sta3Ind,2),'.');
ylabel('Range-Rate [km/s]')
xlabel('Time (hours)')

ind_lastOb = find(Xtrue==options.end_time);
state_errors = X_bestfit(:,1:6) - Xtrue(1:ind_lastOb,2:7);
state_time = [options.start_time:options.inc_time:options.end_time];

figure();
set(gcf,'Color','w');
subplot(3,1,1);
ii = 1;
plot(state_time./(3600),state_errors(:,ii),'.',state_time./(3600),3.*sqrt(Phist(:,(ii-1)*6 + ii)),'r',state_time./(3600),-3.*sqrt(Phist(:,(ii-1)*6 + ii)),'r');
title('ECI Position Errors')
ylabel('X-ECI Errors [km]');
subplot(3,1,2);
ii = 2;
plot(state_time./(3600),state_errors(:,ii),'.',state_time./(3600),3.*sqrt(Phist(:,(ii-1)*6 + ii)),'r',state_time./(3600),-3.*sqrt(Phist(:,(ii-1)*6 + ii)),'r');
ylabel('Y-ECI Errors [km]');
subplot(3,1,3);
ii = 3;
plot(state_time./(3600),state_errors(:,ii),'.',state_time./(3600),3.*sqrt(Phist(:,(ii-1)*6 + ii)),'r',state_time./(3600),-3.*sqrt(Phist(:,(ii-1)*6 + ii)),'r');
ylabel('Z-ECI Errors [km]');
xlabel('Time (hours)')

figure();
set(gcf,'Color','w');
subplot(3,1,1);
ii = 4; 
plot(state_time./(3600),state_errors(:,ii),'.',state_time./(3600),3.*sqrt(Phist(:,(ii-1)*6 + ii)),'r',state_time./(3600),-3.*sqrt(Phist(:,(ii-1)*6 + ii)),'r');
title('ECI Velocity Errors')
ylabel('Xd-ECI Errors [km/s]');
subplot(3,1,2);
ii = 5;
plot(state_time./(3600),state_errors(:,ii),'.',state_time./(3600),3.*sqrt(Phist(:,(ii-1)*6 + ii)),'r',state_time./(3600),-3.*sqrt(Phist(:,(ii-1)*6 + ii)),'r');
ylabel('Yd-ECI Errors [km/s]');
subplot(3,1,3);
ii = 6;
plot(state_time./(3600),state_errors(:,ii),'.',state_time./(3600),3.*sqrt(Phist(:,(ii-1)*6 + ii)),'r',state_time./(3600),-3.*sqrt(Phist(:,(ii-1)*6 + ii)),'r');
ylabel('Zd-ECI Errors [km/s]');
xlabel('Time (hours)')

pos_err_RSW = zeros(length(state_errors),3);
PDiag_RSW = zeros(length(state_errors),3);
for i = 1:length(state_errors)
    Rvec = Xtrue(i,2:4)./norm(Xtrue(i,2:4));
    Wvec = cross(Xtrue(i,2:4),Xtrue(i,5:7))./norm(cross(Xtrue(i,2:4),Xtrue(i,5:7)));
    Svec = cross(Wvec,Rvec);
    ECI2RSW = [Rvec; Svec; Wvec];
    pos_err_RSW(i,:) = (ECI2RSW*state_errors(i,1:3)')';
    tempP = reshape(Phist(i,:),6,6);
    PDiag_RSW(i,1:3) = diag(ECI2RSW*tempP(1:3,1:3)*ECI2RSW');
%     PDiag_RSW(i,4:6) = (ECI2RSW*Phist(i,4:6)')';
end

figure();
set(gcf,'Color','w');
subplot(3,1,1);
plot(state_time./(3600),pos_err_RSW(:,1),'.',state_time./(3600),3.*sqrt(PDiag_RSW(:,1)),'r',state_time./(3600),-3.*sqrt(PDiag_RSW(:,1)),'r');
title('RSW Position Errors')
ylabel('X-RSW Errors [km]');
subplot(3,1,2);
plot(state_time./(3600),pos_err_RSW(:,2),'.',state_time./(3600),3.*sqrt(PDiag_RSW(:,2)),'r',state_time./(3600),-3.*sqrt(PDiag_RSW(:,2)),'r');
ylabel('Y-RSW Errors [km]');
subplot(3,1,3);
plot(state_time./(3600),pos_err_RSW(:,3),'.',state_time./(3600),3.*sqrt(PDiag_RSW(:,3)),'r',state_time./(3600),-3.*sqrt(PDiag_RSW(:,3)),'r');
ylabel('Z-RSW Errors [km]');
xlabel('Time (hours)')

RMS_x = sqrt(sum(state_errors(RMS_zero_ind:end,1).*state_errors(RMS_zero_ind:end,1))/length(state_errors(RMS_zero_ind:end,1)));
RMS_y = sqrt(sum(state_errors(RMS_zero_ind:end,2).*state_errors(RMS_zero_ind:end,2))/length(state_errors(RMS_zero_ind:end,1)));
RMS_z = sqrt(sum(state_errors(RMS_zero_ind:end,3).*state_errors(RMS_zero_ind:end,3))/length(state_errors(RMS_zero_ind:end,1)));
RMS_xdot = sqrt(sum(state_errors(RMS_zero_ind:end,4).*state_errors(RMS_zero_ind:end,4))/length(state_errors(RMS_zero_ind:end,1)));
RMS_ydot = sqrt(sum(state_errors(RMS_zero_ind:end,5).*state_errors(RMS_zero_ind:end,5))/length(state_errors(RMS_zero_ind:end,1)));
RMS_zdot = sqrt(sum(state_errors(RMS_zero_ind:end,6).*state_errors(RMS_zero_ind:end,6))/length(state_errors(RMS_zero_ind:end,1)));

RMS_xRSW = sqrt(sum(pos_err_RSW(RMS_zero_ind:end,1).*pos_err_RSW(RMS_zero_ind:end,1))/length(pos_err_RSW(RMS_zero_ind:end,1)));
RMS_yRSW = sqrt(sum(pos_err_RSW(RMS_zero_ind:end,2).*pos_err_RSW(RMS_zero_ind:end,2))/length(pos_err_RSW(RMS_zero_ind:end,1)));
RMS_zRSW = sqrt(sum(pos_err_RSW(RMS_zero_ind:end,3).*pos_err_RSW(RMS_zero_ind:end,3))/length(pos_err_RSW(RMS_zero_ind:end,1)));

fprintf('\n---------------------------------------------\n\n');
fprintf('POS RMS  = %g %g %g km\n', RMS_x, RMS_y, RMS_z);
fprintf('VEL RMS  = %g %g %g km/s\n', RMS_xdot, RMS_ydot, RMS_zdot);
fprintf('POS RSW RMS  = %g %g %g km\n', RMS_xRSW, RMS_yRSW, RMS_zRSW);
fprintf('\n---------------------------------------------\n\n');

if useFilter > 1
    
    state_errors_est = zeros(length(X_pf),6);
    for ii = 1:length(obs.times)
        useInd = find(Xtrue(:,1)==obs.times(ii),1);
        state_errors_est(ii,:) = X_pf(ii,1:6) - Xtrue(useInd,2:7);
    end
    
    figure();
    set(gcf,'Color','w');
    subplot(3,1,1);
    ii = 1;
    plot(obs.times./(3600),state_errors_est(:,ii),'.',obs.times./(3600),3.*sqrt(P_pf(:,(ii-1)*6 + ii)),'r',obs.times./(3600),-3.*sqrt(P_pf(:,(ii-1)*6 + ii)),'r');
    title('KF xhat ECI Position Errors')
    ylabel('X-ECI Errors [km]');
    subplot(3,1,2);
    ii = 2;
    plot(obs.times./(3600),state_errors_est(:,ii),'.',obs.times./(3600),3.*sqrt(P_pf(:,(ii-1)*6 + ii)),'r',obs.times./(3600),-3.*sqrt(P_pf(:,(ii-1)*6 + ii)),'r');
    ylabel('Y-ECI Errors [km]');
    subplot(3,1,3);
    ii = 3;
    plot(obs.times./(3600),state_errors_est(:,ii),'.',obs.times./(3600),3.*sqrt(P_pf(:,(ii-1)*6 + ii)),'r',obs.times./(3600),-3.*sqrt(P_pf(:,(ii-1)*6 + ii)),'r');
    ylabel('Z-ECI Errors [km]');
    xlabel('Time (hours)')
    
    figure();
    set(gcf,'Color','w');
    subplot(3,1,1);
    ii = 4;
    plot(obs.times./(3600),state_errors_est(:,ii),'.',obs.times./(3600),3.*sqrt(P_pf(:,(ii-1)*6 + ii)),'r',obs.times./(3600),-3.*sqrt(P_pf(:,(ii-1)*6 + ii)),'r');
    title('KF xhat ECI Velocity Errors')
    ylabel('Xd-ECI Errors [km/s]');
    subplot(3,1,2);
    ii = 5;
    plot(obs.times./(3600),state_errors_est(:,ii),'.',obs.times./(3600),3.*sqrt(P_pf(:,(ii-1)*6 + ii)),'r',obs.times./(3600),-3.*sqrt(P_pf(:,(ii-1)*6 + ii)),'r');
    ylabel('Yd-ECI Errors [km/s]');
    subplot(3,1,3);
    ii = 6;
    plot(obs.times./(3600),state_errors_est(:,ii),'.',obs.times./(3600),3.*sqrt(P_pf(:,(ii-1)*6 + ii)),'r',obs.times./(3600),-3.*sqrt(P_pf(:,(ii-1)*6 + ii)),'r');
    ylabel('Zd-ECI Errors [km/s]');
    xlabel('Time (hours)')
    
    % Plot diagonals of covariance on log plot
    figure; 
    set(gcf,'Color','w');
    for ii = 1:6
        semilogy(obs.times./(3600),sqrt(P_pf(:,(ii-1)*6 + ii)),'-o');
        hold on;
    end
    ylabel('\sigma [km, km/s]');
    xlabel('Time (hours)');
    legend('x','y','z','xd','yd','zd');
    
    % Plot trace of covariance
    figure; 
    set(gcf,'Color','w');
    inds = zeros(1,6);
    for ii = 1:6
        inds(ii) = (ii-1)*6 + ii;
    end
    semilogy(obs.times./3600, sum(P_pf(:,inds),2),'-o')
    ylabel('Trace of P');
    xlabel('Time (hours)');

end

toc