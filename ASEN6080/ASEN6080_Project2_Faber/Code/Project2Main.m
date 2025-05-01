%% ASEN 6080 Project 2 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
    % Path logistics
addpath("..\")
addpath(genpath("..\..\Utilities\"));

    % Celestial body constants
pConst = getPlanetConst();

    % Spacecraft constants
scConst = getSCConst(2);

    % Initialize stations struct
stations = makeStations(pConst, 2);

    % ode45 options
opt = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

    %% Ask user for input
menuText = sprintf("\n\t\tWelcome to Ian Faber's ASEN 6080 Project 2 Code!" + ...
                   "\nChoose an option to do by entering the corresponding number:" + ...
                   "\n\t1. Run Part 1: Dynamics, Measurements, and B Plane Verification (estimated run time: ~10 seconds)" + ...
                   "\n\t2. Run Part 2: Estimating State with Known Target and Models (estimated run time: ~1.5 minutes)" + ...
                   "\n\t3. Run Part 3: Fitting Data with Unknown Issues and No Truth (estimated run time: ~7 minutes)" + ...
                   "\n\t4. Exit (estimated run time: instant ;-) )" + ...
                   "\nChoice: ");
choice = input(menuText);

switch choice
    case 1
        %% Part 1: Dynamics Implementation
        fprintf("\nProving Dynamics, Measurement, and B Plane Implementations\n")
        
            % Set up initial state and get initial accelerations
        X0 = [scConst.X0_cart; scConst.C_R];
        aTest = orbitEOM_MuSunSRP(0, X0, pConst, scConst);
        
            % Load truth trajectory to verify dynamics are working
        truthTraj = load("..\Data\Project2_Prob2_truth_traj_50days.mat");
        
            % Integrate X0 to verify EOM function
        tspan = [truthTraj.Tt_50; (truthTraj.Tt_50(end)+1000:1000:200*24*60*60)']; % Add time up to 200 days in sec
        [t_test, X_test] = ode45(@(t,X)orbitEOM_MuSunSRP(t,X,pConst,scConst), tspan, truthTraj.Xt_50(1,1:7)', opt);
        
            % Plot error between modeled and true trajectory
        idx = 1:length(truthTraj.Tt_50); % Only plot data when truth data exists
        
        figure; tl = tiledlayout(2,3);
        title(tl, "Part 1: Error between modeled and true trajectory")
        nexttile;
            hold on; grid on;
            plot(t_test(idx), X_test(idx,1) - truthTraj.Xt_50(idx,1), 'b-');
            xlabel("Time [sec]"); ylabel("\DeltaX [km]")
        nexttile;
            hold on; grid on;
            plot(t_test(idx), X_test(idx,2) - truthTraj.Xt_50(idx,2), 'b-');
            xlabel("Time [sec]"); ylabel("\DeltaY [km]")
        nexttile;
            hold on; grid on;
            plot(t_test(idx), X_test(idx,3) - truthTraj.Xt_50(idx,3), 'b-');
            xlabel("Time [sec]"); ylabel("\DeltaZ [km]")
        nexttile;
            hold on; grid on;
            plot(t_test(idx), X_test(idx,4) - truthTraj.Xt_50(idx,4), 'b-');
            xlabel("Time [sec]"); ylabel("\DeltaXdot [km/s]")
        nexttile;
            hold on; grid on;
            plot(t_test(idx), X_test(idx,5) - truthTraj.Xt_50(idx,5), 'b-');
            xlabel("Time [sec]"); ylabel("\DeltaYdot [km/s]")
        nexttile;
            hold on; grid on;
            plot(t_test(idx), X_test(idx,6) - truthTraj.Xt_50(idx,6), 'b-');
            xlabel("Time [sec]"); ylabel("\DeltaZdot [km/s]")
        
            % Plot initial situation
        [R_Earth, ~, ~, ~] = Ephem(pConst.initEpoch, 3, 'EME2000');
        R_E_t = []; R_E_t_2 = [];
        for k = 1:length(t_test)%truthTraj.Tt_50)
            % [R_E_t(:,k), ~, ~, ~] = Ephem(pConst.initEpoch + truthTraj.Tt_50(k)/(24*60*60), 3, 'EME2000'); % Earth orbit over this scenario
            [R_E_t(:,k), ~, ~, ~] = Ephem(pConst.initEpoch + t_test(k)/(24*60*60), 3, 'EME2000'); % Earth orbit over this scenario
        end
        
        for k = 1:365
            [R_E_t_2(:,k), ~, ~, ~] = Ephem(pConst.initEpoch + k, 3, 'EME2000'); % Earth full orbit
        end
        
        figure;
        title("Part 1: Truth orbit comparison until 200 days - Sun origin frame (EME2000)")
        hold on; grid on; axis equal;
        plot3(0,0,0,'y.','MarkerSize', 50);
        plot3(X_test(1,1)+R_Earth(1), X_test(1,2)+R_Earth(2), X_test(1,3)+R_Earth(3),'m.','MarkerSize',30);
        plot3(truthTraj.Xt_50(:,1)+R_E_t(1,1:length(truthTraj.Tt_50))', truthTraj.Xt_50(:,2)+R_E_t(2,1:length(truthTraj.Tt_50))', truthTraj.Xt_50(:,3)+R_E_t(3,1:length(truthTraj.Tt_50))', 'm-.', 'LineWidth', 3)
        plot3(X_test(:,1)+R_E_t(1,:)', X_test(:,2)+R_E_t(2,:)', X_test(:,3)+R_E_t(3,:)', 'm:', 'LineWidth', 2)
        plot3(R_Earth(1), R_Earth(2), R_Earth(3),'b.','MarkerSize',40);
        plot3(R_E_t(1,:), R_E_t(2,:), R_E_t(3,:), 'b-.', 'LineWidth', 2)
        plot3(R_E_t_2(1,:), R_E_t_2(2,:), R_E_t_2(3,:), 'k:')
        xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]")
        legend("Sun", "SC", "SC Truth Orbit", "SC Modeled Orbit", "Earth", "Earth Orbit")
        view([30 35])
        
        figure;
        title("Part 1: Truth orbit comparison until 200 days - Earth origin frame (ECI)")
        hold on; grid on; axis equal;
        plot3(0,0,0,'b.','MarkerSize', 50);
        plot3(X_test(1,1), X_test(1,2), X_test(1,3),'m.','MarkerSize',30);
        plot3(truthTraj.Xt_50(:,1), truthTraj.Xt_50(:,2), truthTraj.Xt_50(:,3), 'm-.', 'LineWidth', 3)
        plot3(X_test(:,1), X_test(:,2), X_test(:,3), 'm:', 'LineWidth', 2)
        xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]")
        legend("Earth", "SC", "SC Truth Orbit", "SC Modeled Orbit")
        view([30 35])
        
        %% Part 1: DSN Measurement Processing
        data2a = readmatrix("..\Data\Project2a_Obs.txt");
        [tMeas, stations] = readObsData(stations, data2a);
        
        tspan = tMeas; % Add time up to 200 days in sec
        [t_test, X_test] = ode45(@(t,X)orbitEOM_MuSunSRP(t,X,pConst,scConst), tspan, truthTraj.Xt_50(1,1:7)', opt);
        
            % Propagate station states and generate clean measurements
        for k = 1:length(tMeas)
            dTheta = pConst.wEarth*tMeas(k);
            for kk = 1:length(stations)
                r = rotZ(dTheta)*stations(kk).X0;
                v = cross([0;0;pConst.wEarth], r);
        
                stations(kk).Xs = [stations(kk).Xs; [r', v', tMeas(k)]];        
            end
        end
        
        stations_nom = makeStations(pConst, 2);
            % Propagate station states and generate clean measurements
        for k = 1:length(t_test)
            dTheta = pConst.wEarth*t_test(k);
            for kk = 1:length(stations_nom)
                r = rotZ(dTheta)*stations_nom(kk).X0;
                v = cross([0;0;pConst.wEarth], r);
        
                stations_nom(kk).Xs = [stations_nom(kk).Xs; [r', v', t_test(k)]];
                
                y = generateRngRngRate(X_test(k,:), stations_nom(kk).Xs(k,:), stations_nom(kk).elMask);
        
                if ~isnan(y)
                    stations_nom(kk).rho = [stations_nom(kk).rho; y(1)];
                    stations_nom(kk).rhoDot = [stations_nom(kk).rhoDot; y(2)];
                    stations_nom(kk).elAngle = [stations_nom(kk).elAngle; y(3)];
                    stations_nom(kk).tMeas = [stations_nom(kk).tMeas; t_test(k)];
                    if true
                        R = diag([stations_nom(kk).sigRho^2, stations_nom(kk).sigRhoDot^2]);
                        stations_nom(kk).R = [stations_nom(kk).R; {R}];
                    else
                        R = zeros(2,2);
                        stations_nom(kk).R = [stations_nom(kk).R; {R}];
                    end
                end
            end
        end
        
        titleText = sprintf("Part 1: Provided Measurement Data");
        xLabel = "Time [sec]"; 
        yLabel = ["\rho [km]", "\rhoDot [km/s]"];
        plotMeasurements(stations, titleText, xLabel, yLabel);
        
        titleText = sprintf("Part 1: Propagated Measurement Data");
        xLabel = "Time [sec]"; 
        yLabel = ["\rho [km]", "\rhoDot [km/s]"];%, "Elevation Angle [deg]"];
        plotMeasurements(stations_nom, titleText, xLabel, yLabel);
        
        %% Part 1: Bplane implementation
            % Set up ODE events for sphere of influence
        opt.Events = @(t,X)SOICheck(t,X,pConst);
        
            % Define function handle for integration
        DynFunc = @(t,XPhi)STMEOM_MuSunSRP(t,XPhi,pConst,scConst);
        
            % Integrate to 3 SOI
        XPhi_0 = [truthTraj.Xt_50(1,1:7)'; reshape(eye(7), 49, 1)];
        P0 = 1e-11*eye(7); % Test covariance - not physical
        tspan_3SOI = [truthTraj.Tt_50; (truthTraj.Tt_50(end)+1000:1000:300*24*60*60)']; % Add time up to 300 days in sec
        [t_3SOI, XPhi_3SOI] = ode45(@(t,XPhi)DynFunc(t,XPhi), tspan_3SOI, XPhi_0, opt);
        Phi_3SOI = reshape(XPhi_3SOI(end, 8:end),7,7);
        P_3SOI = Phi_3SOI*P0*Phi_3SOI';
        
            % Calculate Bplane without SOI checker
        [BdotR_truth, BdotT_truth, ~, ~, ~, X_crossing_truth, P_BPlane, STR2ECI, XPhi_BPlane, t_BPlane] = calcBPlane(XPhi_3SOI(end,:)', t_3SOI(end), P_3SOI, pConst, DynFunc, odeset('RelTol',1e-13,'AbsTol',1e-13));
        
            % Plot BdotR and BdotT location + uncertainty ellipse
        boundLevel = 3;
        titleText = sprintf("Part 1: B Plane Target Estimate from Truth Data");
        xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; 
        ellipseLabel = sprintf("\\pm %.0f\\sigma B plane target uncertainty, dummy P_0", boundLevel);
        ellipseColor = 'b';
        BVecLabel = sprintf("B Plane Target after propagating truth data: \nBdotR = %.4e km,\nBdotT = %.4e km", BdotR_truth, BdotT_truth);
        plotBPlane(BdotR_truth, BdotT_truth, X_crossing_truth, P_BPlane, STR2ECI, pConst, boundLevel, titleText, xLabel, yLabel, zLabel, ellipseLabel, ellipseColor, BVecLabel, 420, true, 1);
        
        figure(420)
        nexttile(1)
            hold on;
            idx = length(t_BPlane); offset = 1000;
            orbit = plot3(XPhi_BPlane(idx-offset:end,1), XPhi_BPlane(idx-offset:end,2), XPhi_BPlane(idx-offset:end,3), 'm--', 'DisplayName', "Modeled Orbit");
        
    case 2   
        %% Part 2: Estimate State with Known Target and Models with varying amounts of data
        fprintf("\nEstimating State with Known Target and Models\n")
        
            % Set initial state estimate and deviation
        X0 = [scConst.X0_cart; scConst.C_R];
        % X0 = truthTraj.Xt_50(end,1:7)';
        x0 = zeros(size(X0));

            % Set initial covariances
        sigR = 100; % km
        sigV = 0.1; % km/s
        sigCR = 0.1;
        P0 = diag([sigR^2, sigR^2, sigR^2, sigV^2, sigV^2, sigV^2, sigCR^2]);
        
            % Reset ODE45 options - won't need SOI checker
        opt = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
        
            % Load truth trajectory
        truthTraj = load("..\Data\Project2_Prob2_truth_traj_50days.mat");

            % Remake and populate stations structure
                % Make structure
        stations = makeStations(pConst, 2);
                % Read and populate observations
        data2a = readmatrix('..\Data\Project2a_Obs.txt');
        [tMeas, stations] = readObsData(stations, data2a);
        
                % Propagate station states
        for k = 1:length(tMeas)
            dTheta = pConst.wEarth*tMeas(k);
            for kk = 1:length(stations)
                r = rotZ(dTheta)*stations(kk).X0;
                v = cross([0;0;pConst.wEarth], r);
        
                stations(kk).Xs = [stations(kk).Xs; [r', v', tMeas(k)]];        
            end
        end
        
            % Plot measurements
        titleText = sprintf("Part 2: Provided Measurement Data");
        xLabel = "Time [sec]"; 
        yLabel = ["\rho [km]", "\rhoDot [km/s]"];
        plotMeasurements(stations, titleText, xLabel, yLabel);
        
            % Loop over days of interest
        days = [50; 100; 150; 200];
        dayColors = ['b', 'r', 'c', 'k'];
        LKFRuns = [];
        EKFRuns = [];
        t_Combined = [];
        X_Combined = [];
        PEst_Combined = {};
        prefits_Combined = [];
        postfits_Combined = [];
        
        XPhi_0 = [truthTraj.Xt_50(1,1:7)'; reshape(eye(length(X0)),length(X0)^2,1)];
        [t_days, XPhi_days] = ode45(@(t,X)STMEOM_MuSunSRP(t,X,pConst,scConst), tMeas, XPhi_0, opt);
        X_days = XPhi_days(:,1:7);
        Phi_days = XPhi_days(:,8:end);
        Phi_ext = {};
        for k = 1:length(t_days)
            Phi_ext = [Phi_ext; reshape(Phi_days(k,:), 7, 7)];
        end

        for k = 1:length(days)
            fprintf("\n--- Filtering from %.0f to %.0f days of data ---\n", days(k) - 50, days(k));
        
                 % Properly initialize filters
            if k == 1
                X0_i = X0;
                x0_i = x0;
                P0_i = P0;
            else
                X0_i = X_days(kEnd_truth + 1, :)';
                x0_i = x0;

                     % Find STM to propagate X and P between sections
                        % Filter section
                Phi_tStartEnd = LKFRuns(end).LKFOut.Phi_full;
                Phi_tStartEnd = Phi_tStartEnd(1:7,1:7);
                        % Prefilter section
                Phi_t0tStart = Phi_ext{kStart};
                        % Overall time STM from t0 to t
                Phi_t0t = Phi_ext{kEnd_truth + 1};
                Phi = (Phi_tStartEnd^-1)*(Phi_t0tStart^-1)*Phi_t0t;
                    % Propagate P
                P = PEst_Combined{end};
                P = P(1:7,1:7);
                P0_i = 10000*Phi*P*Phi'; % Inflate propagated uncertainty a bit to avoid saturation

                % P0_i = 0.5*P0; % Start new segment with high uncertainty to avoid saturation
            end


                % Create new truth data for the specified number of days
            if k == 1
                kStart = 1;   
            else
                kStart = kEnd_truth + 1;
            end
            kEnd_truth = find(tMeas <= days(k)*24*60*60, 1, 'last'); 
        
            % XPhi_0 = [truthTraj.Xt_50(1,1:7)'; reshape(eye(length(X0)),length(X0)^2,1)];
            % [t_days, XPhi_days] = ode45(@(t,X)STMEOM_MuSunSRP(t,X,pConst,scConst), tMeas(1:kEnd_truth), XPhi_0, opt);
            % X_days = XPhi_days(:,1:7);
        
           
        
            %     % Determine cutoffs between LKF and EKF, i.e. when a large measurement
            %     % gap happens
            % dt = diff(t_days); % Find time between measurements
            % med_dt = median(dt); % Find median time gap
            % kSec = [find(dt > 1000*med_dt); kEnd_truth]; % Pull out when time gaps are larger than the median
            
                % Configure plotting
                    % Only plot residuals and state errors w/ bounds
            plotBool = [false; false; false; false; false; false; false; false];     
            
                % Run LKF on data
            numMeas = kEnd_truth - kStart;
            LKFRun = runLKF(X0_i, x0_i, P0_i, pConst, scConst, stations, X_days(kStart:kEnd_truth,:), t_days(kStart:kEnd_truth), tMeas(kStart), numMeas, 10, plotBool);
            LKFRuns = [LKFRuns; LKFRun];
        
                % Accumulate data
            t_Combined = [t_Combined; LKFRun.t_LKF];
            X_Combined = [X_Combined, LKFRun.X_LKF];
            PEst_Combined = [PEst_Combined, LKFRun.LKFOut.PEst];
            prefits_Combined = [prefits_Combined, LKFRun.LKFOut.prefit_res];
            postfits_Combined = [postfits_Combined, LKFRun.LKFOut.postfit_res];
        
            %     % Run LKF/EKF, initializing with LKF after a large gap
            % n = 100000; 
            % if k == 1
            %     start = 1;
            % else
            %     start = sec + 1;
            % end
            % for sec = start:length(kSec)
            %     kEnd = kSec(sec);
            % 
            %     if tMeas(kStart) == tMeas(kEnd) % Only 1 measurement
            %         continue;
            %     end
            % 
            %     fprintf("\n---- Filtering from t = %.3f to t  = %.3f ----\n", tMeas(kStart), tMeas(kEnd));
            % 
            %         % Determine when to turn on EKF
            %     if kEnd - kStart < n % Use only LKF
            %         numMeas = kEnd - kStart;
            %     else % Initialize with n LKF measurements, then switch to EKF
            %         numMeas = n;
            %     end
            % 
            %         % Run LKF
            %     LKFRun = runLKF(X0_i, x0_i, P0_i, pConst, scConst, stations, X_days(kStart:kEnd,:), t_days(kStart:kEnd), tMeas(kStart), numMeas, 10, plotBool);
            %     LKFRuns = [LKFRuns; LKFRun];
            %     P0_i = LKFRun.LKFOut.PEst{end};
            % 
            %         % Accumulate data
            %     t_Combined = [t_Combined; LKFRun.t_LKF];
            %     X_Combined = [X_Combined, LKFRun.X_LKF];
            %     PEst_Combined = [PEst_Combined, LKFRun.LKFOut.PEst];
            %     prefits_Combined = [prefits_Combined, LKFRun.LKFOut.prefit_res];
            %     postfits_Combined = [postfits_Combined, LKFRun.LKFOut.postfit_res];
            % 
            %         % Run LKF if we have enough measurments to start
            %     if numMeas == n
            %         t_start = LKFRun.t_LKF(end); t_end = tMeas(kEnd);
            % 
            %         if t_start == t_end % Only 1 measurement
            %             continue;
            %         end
            % 
            %         X_start = LKFRun.X_LKF(:,end) - LKFRun.LKFOut.xEst(:,end);
            %         P_start = LKFRun.LKFOut.PEst{end}; 
            % 
            %         EKFRun = runEKF(X_start, P_start, pConst, scConst, stations, X_days(kStart:kEnd,:), t_days(kStart:kEnd), t_start, t_end, plotBool);
            %         EKFRuns = [EKFRuns; EKFRun];
            % 
            %         P0_i = EKFRun.EKFOut.PEst{end};
            % 
            %                 % Accumulate data
            %         t_Combined = [t_Combined; EKFRun.t_EKF];
            %         X_Combined = [X_Combined, EKFRun.X_EKF];
            %         PEst_Combined = [PEst_Combined, EKFRun.EKFOut.PEst];
            %         prefits_Combined = [prefits_Combined, EKFRun.EKFOut.prefit_res];
            %         postfits_Combined = [postfits_Combined, EKFRun.EKFOut.postfit_res];
            %     end
            % 
            %     %     % Run IEKF
            %     % EKFRun = runIEKF(X0_i,P0_i,pConst,scConst,stations,X_days(kStart:kEnd,:),t_days(kStart:kEnd),tMeas(kStart),tMeas(kEnd),plotBool);
            %     % EKFRuns = [EKFRuns; EKFRun];
            % 
            %     %     % Accumulate data
            %     % t_Combined = [t_Combined; EKFRun.t_EKF];
            %     % X_Combined = [X_Combined, EKFRun.X_EKF];
            %     % PEst_Combined = [PEst_Combined, EKFRun.EKFOut.PEst];
            %     % prefits_Combined = [prefits_Combined, EKFRun.EKFOut.prefit_res];
            %     % postfits_Combined = [postfits_Combined, EKFRun.EKFOut.postfit_res];
            % 
            %         % Reset for next set of measurements
            %     if sec < length(kSec)
            %             % Update starting index
            %         kStart = kEnd + 1;
            % 
            %             % Get starting covariance, Xstar, and deviation for next time
            %         Phi = eye(length(X0));
            %         XPhi_0 = [X_Combined(:,end); reshape(Phi, length(X0)^2,1)];
            %         [~,XPhi] = ode45(@(t,XPhi)STMEOM_MuSunSRP(t,XPhi,pConst,scConst), [t_Combined(end), t_days(kStart)], XPhi_0, opt);
            %         Phi = reshape(XPhi(end,(length(X0)+1):end), length(X0), length(X0));
            %         X = XPhi(end,1:length(X0))';
            % 
            %         P0_i = Phi*PEst_Combined{end}*Phi';
            % 
            %         % X0_i = X_days(kStart,:)';
            %         % X0_i = X;
            %         % x0_i = zeros(size(X0));
            % 
            %         %     % Update X0 and x0
            %         % if numMeas == n % EKF was used, Xstar changed
            %         %     [~,X_prop] = ode45(@(t,X)orbitEOM_MuSunSRP(t,X,pConst,scConst), [t_Combined(end), t_days(kStart)], X_Combined(:,end), opt);
            %         %     X0_EKF = X_prop(end,:)';
            %         %     X0_i = X0_EKF;
            %         %     x0_i = zeros(size(X0)); %X0_EKF - X_days(kStart,:)';
            %         % else % Only CKF was used
            %         %     X0_i = X_days(kStart,:)';
            %         %     x0_i = zeros(size(X0));
            %         % end
            % 
            %         X0_i = X_days(kStart,:)';
            %         x0_i = zeros(size(X0));
            %     end
            % end
            
                % Calculate new X_ref vector
            X_ref_Combined = [];
            for kk = 1:length(t_Combined)
                X_ref_Combined = [X_ref_Combined; X_days(t_days == t_Combined(kk),:)];
            end
            
                % Calculate state error and uncertainty
            stateError_Combined = X_Combined' - X_ref_Combined;
            
            sigma_Combined = [];
            for kk = 1:length(PEst_Combined)
                P = PEst_Combined{kk};
            
                sigPart = [];
                for elem = 1:size(P,1)
                    sigPart = [sigPart, sqrt(P(elem,elem))];
                end
            
                sigma_Combined = [sigma_Combined; sigPart];
            end
            
            titleText = sprintf("Pre-Fit Residuals over %.0f days", days(k)); 
            xLabel = "Time [sec]"; 
            yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
            colors = ['b', 'r'];
            plotResiduals(t_Combined, prefits_Combined, titleText, xLabel, yLabel, colors);
        
            titleText = sprintf("Post-Fit Residuals over %.0f days", days(k)); 
            xLabel = "Time [sec]"; 
            yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
            colors = ['b', 'r'];
            plotResiduals(t_Combined, postfits_Combined, titleText, xLabel, yLabel, colors);
        
            titleText = sprintf("Part 2: State error over %.0f days", days(k));
            xLabel = "Time [sec]";
            yLabel = ["X error [km]", "Y error [km]", "Z error [km]", ...
                      "Xdot error [km/s]", "Ydot error [km/s]", "Zdot error [km/s]", ...
                      "C_R error [n.d.]"];
            plotStateError(t_Combined, stateError_Combined, t_Combined, sigma_Combined, 3, titleText, xLabel, yLabel);
            
            % titleText = sprintf("Part 2: Nominal vs. Estimated State over %.3f days", t_Combined(end)/(24*60*60));
            % xLabel = "Time [sec]";
            % yLabel = ["X [km]", "Y [km]", "Z [km]", ...
            %           "Xdot [km/s]", "Ydot [km/s]", "Zdot [km/s]", ...
            %           "C_R [n.d.]"];
            % plotStateError(t_Combined, X_Combined', [], [], [], titleText, xLabel, yLabel);
            % for kk = 1:size(X_days,2)
            %     nt = nexttile(kk);
            %         hold on; grid on;
            %         plot(t_days, X_days(:,kk),'k--','DisplayName',"Nominal Orbit");
            % end
            
            %% Part 2: Bplane Implementation
            fprintf("\n---- Calculating B Plane for %.0f days ----\n", days(k))
            
                % Set up ODE events for sphere of influence
            opt.Events = @(t,X)SOICheck(t,X,pConst);
            
                % Define function handle for integration
            DynFunc = @(t,XPhi)STMEOM_MuSunSRP(t,XPhi,pConst,scConst);
            
                % Integrate to 3 SOI
            XPhi_0 = [X_Combined(:,end); reshape(eye(7), 49, 1)];
            P0_Bplane = PEst_Combined{end};
            tspan_3SOI = (t_Combined(end):1000:300*24*60*60)'; % Add time up to 300 days in sec
            [t_3SOI, XPhi_3SOI] = ode45(@(t,XPhi)DynFunc(t,XPhi), tspan_3SOI, XPhi_0, opt);
            Phi_3SOI = reshape(XPhi_3SOI(end, 8:end),7,7);
            P_3SOI = Phi_3SOI*P0_Bplane*Phi_3SOI';
            
                % Calculate Bplane without SOI checker
            [BdotR, BdotT, ~, ~, ~, X_crossing, P_BPlane, STR2ECI, XPhi_BPlane, t_BPlane] = calcBPlane(XPhi_3SOI(end,:)', t_3SOI(end), P_3SOI, pConst, DynFunc, odeset('RelTol',1e-13,'AbsTol',1e-13));
            
                % Plot BdotR and BdotT location + uncertainty ellipse
            if k == 1
                newFig = true;
            else
                newFig = false;
            end
        
            boundLevel = 3;
            titleText = sprintf("Part 2: B Plane Target Estimate from Estimated State");
            xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; 
            ellipseLabel = sprintf("%.0f\\sigma B plane target uncertainty, %.3f days of data", boundLevel, t_Combined(end)/(24*60*60));
            ellipseColor = dayColors(k);
            BVecLabel = sprintf("\nB Plane Target after %.3f days of data: \nBdotR = %.4e km,\nBdotT = %.4e km", t_Combined(end)/(24*60*60), BdotR, BdotT);
            plotBPlane(BdotR, BdotT, X_crossing, P_BPlane, STR2ECI, pConst, boundLevel, titleText, xLabel, yLabel, zLabel, ellipseLabel, ellipseColor, BVecLabel, 69, newFig, 2);
        
        end

            % Truth target
        BdotR_truth = 14970.824; 
        BdotT_truth = 9796.737;
        xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; 
        ellipseLabel = sprintf("Scaled %.0f\\sigma B plane target uncertainty from 200 days of data", boundLevel);
        ellipseColor = 'g';
        BVecLabel = sprintf("\nB Plane Target from truth data: \nBdotR = %.4e km,\nBdotT = %.4e km", BdotR_truth, BdotT_truth);
        plotBPlane(BdotR_truth, BdotT_truth, X_crossing, 1e-9*P_BPlane, STR2ECI, pConst, boundLevel, titleText, xLabel, yLabel, zLabel, ellipseLabel, ellipseColor, BVecLabel, 69, false, 2);
    case 3
        %% Part 3: Fit Data with Unknown Issues and No Truth
        fprintf("\nFitting Data with Unknown Issues and no Truth\n")

            % Update scConst for part 3
        scConst = getSCConst(3);
        
            % Set initial state estimate and deviation
        % X0 = [scConst.X0_cart; scConst.C_R; zeros(3,1)];
        X0 = [scConst.X0_cart; scConst.C_R];
        x0 = zeros(size(X0));
        
            % Set initial covariances
        sigR = 100; % km
        sigV = 0.1; % km/s
        sigCR = 0.1;
        P0 = diag([sigR^2, sigR^2, sigR^2, sigV^2, sigV^2, sigV^2, sigCR^2]);
        
            % Reset ODE45 options - won't need SOI checker
        opt = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
        
            % Remake and populate stations structure
        part = 3;
                % Make structure
        stations = makeStations(pConst, part);
                % Read and populate observations
        data2a = readmatrix('..\Data\Project2b_Obs.txt');
        [tMeas, stations] = readObsData(stations, data2a);
        
                % Propagate station states
        for k = 1:length(tMeas)
            dTheta = pConst.wEarth*tMeas(k);
            for kk = 1:length(stations)
                r = rotZ(dTheta)*stations(kk).X0;
                v = cross([0;0;pConst.wEarth], r);
        
                stations(kk).Xs = [stations(kk).Xs; [r', v', tMeas(k)]];        
            end
        end
        
            % Create expected measurements
        tspan = tMeas;
        [t_test, X_test] = ode45(@(t,X)orbitEOM_MuSunSRP(t,X,pConst,scConst), tspan, X0(1:7), opt);
        
        stations_exp = makeStations(pConst, 3);
            % Propagate station states and generate clean measurements
        for k = 1:length(t_test)
            dTheta = pConst.wEarth*t_test(k);
            for kk = 1:length(stations_exp)
                r = rotZ(dTheta)*stations_exp(kk).X0;
                v = cross([0;0;pConst.wEarth], r);
        
                stations_exp(kk).Xs = [stations_exp(kk).Xs; [r', v', t_test(k)]];
                
                y = generateRngRngRate(X_test(k,:), stations_exp(kk).Xs(k,:), stations_exp(kk).elMask);
        
                if ~isnan(y)
                    stations_exp(kk).rho = [stations_exp(kk).rho; y(1)];
                    stations_exp(kk).rhoDot = [stations_exp(kk).rhoDot; y(2)];
                    stations_exp(kk).elAngle = [stations_exp(kk).elAngle; y(3)];
                    stations_exp(kk).tMeas = [stations_exp(kk).tMeas; t_test(k)];
                    if true
                        R = diag([stations_exp(kk).sigRho^2, stations_exp(kk).sigRhoDot^2]);
                        stations_exp(kk).R = [stations_exp(kk).R; {R}];
                    else
                        R = zeros(2,2);
                        stations_exp(kk).R = [stations_exp(kk).R; {R}];
                    end
                end
            end
        end
        
            % Plot measurements
        titleText = sprintf("Part 3: Provided Measurement Data");
        xLabel = "Time [sec]"; 
        yLabel = ["\rho [km]", "\rhoDot [km/s]"];
        fProv = plotMeasurements(stations, titleText, xLabel, yLabel);
        fProv; nexttile(1);
        xlim([0 2.5e7]); ylim([0 3e8]);
        
        titleText = sprintf("Part 3: Expected Measurement Data up to 250 days");
        xLabel = "Time [sec]"; 
        yLabel = ["\rho [km]", "\rhoDot [km/s]"];
        fExp = plotMeasurements(stations_exp, titleText, xLabel, yLabel);
        fExp; nexttile(1); hold on;
        xlim([0 2.5e7]); ylim([0 3e8]);
        
            % Propagate reference orbit
        XPhi_0 = [X0(1:7); reshape(eye(7),7^2,1)];
        [t_ref, XPhi_ref] = ode45(@(t,X)STMEOM_MuSunSRP(t,X,pConst,scConst), tMeas, XPhi_0, opt);
        X_ref = XPhi_ref(:,1:7);
        Phi_ref = XPhi_ref(:,8:end);
        Phi_ext = {};
        for k = 1:length(t_ref)
            Phi_ext = [Phi_ext; reshape(Phi_ref(k,:), 7, 7)];
        end

            % Initialize reference variables
        X_days = X_ref;
        t_days = t_ref;
        % X_comp = X_days';

            % Process data
        % days = linspace(50,250,5); %[50; 100; 150; 200; 250];
        % days = linspace(100,300,3); 
        days = 250;
        dayColors = ['b', 'r', 'c', 'm', 'k'];
        LKFRuns = [];
        EKFRuns = [];
        t_Combined = [];
        X_Combined = [];
        PEst_Combined = {};
        prefits_Combined = [];
        postfits_Combined = [];
        iterations_Combined = [];
        for k = 1:length(days)
            if length(days) == 1
                fprintf("\n--- Filtering from %.0f to %.0f days of data ---\n", days(k) - days(k), days(k));
            else
                fprintf("\n--- Filtering from %.0f to %.0f days of data ---\n", days(k) - median(diff(days)), days(k));
            end

                % Update indices for the current section
            if k == 1
                kStart = 1;   
            else
                kStart = kEnd_truth + 1;
            end
            kEnd_truth = find(tMeas <= days(k)*24*60*60, 1, 'last'); 

                % Properly initialize filters for this section
            if k == 1
                X0_i = X0;
                x0_i = x0;
                P0_i = P0;
            else
                X_last = X_Combined(1:7,end);

                    % Propagate final state estimate forward
                tspan = tMeas((kStart-1):kEnd_truth);
                [t_days, X_days] = ode45(@(t,X)orbitEOM_MuSunSRP(t,X,pConst,scConst), tspan, X_last, opt);

                % X0_i = X_days(kStart, :)';
                X0_i = X_days(2,:)';
                x0_i = zeros(size(X0_i));
                % P0_i = 0.5*P0;
                P = PEst_Combined{end};
                P0_i = 1.5*P(1:7, 1:7);

                %     % Find STM to propagate X and P between sections
                %         % Filter section
                % Phi_tStartEnd = LKFRuns(end).LKFOut.Phi{end};
                % Phi_tStartEnd = Phi_tStartEnd(1:7,1:7);
                %         % Prefilter section
                % Phi_t0tStart = Phi_ext{kStart};
                %         % Overall time STM from t0 to t
                % Phi_t0t = Phi_ext{kEnd_truth + 1};
                % Phi = (Phi_tStartEnd^-1)*(Phi_t0tStart^-1)*Phi_t0t;
                %     % Propagate P
                % P = PEst_Combined{end};
                % P = P(1:7,1:7);
                % P0_i = 1000*Phi*P*Phi';
                % % P0_i = P0; % Start new segment with high uncertainty to avoid saturation
                % 
                %     % Propagate x and make X0
                % x_filt = LKFRuns(end).LKFOut.xEst(1:7,end);
                % x_filt = Phi*x_filt;
                % x_error = LKFRuns(end).X_LKF(1:7,end) - X_days(kEnd_truth,:)';
                % x_error = Phi*x_error;
                % X0_i = X_days(kEnd_truth + 1, :)' + x_error;
                % X0_i = [X0_i; zeros(3,1)];
                
            end
            
            %     % Update indices for the current section
            % if k == 1
            %     kStart = 1;   
            % else
            %     kStart = kEnd_truth + 1;
            % end
            % kEnd_truth = find(tMeas <= days(k)*24*60*60, 1, 'last'); 

                % Configure plotting
                    % Only plot residuals and state errors w/ bounds
                    % Order: [prefits; postfits; R trace; V trace; R
                    %         Ellipsoid; V Ellipsoid; State Error, no
                    %         bounds; State Error, 3 sigma bounds]
            plotBool_LKF = [false; false; false; false; false; false; false; true];     
            plotBool_EKF = [false; false; false; false; false; false; true; false]; 
            
                % Run LKF/IEKF on data
            numMeas = kEnd_truth - kStart;
        
                    % Set up DMC
                        % LKF
            tau_x_LKF = 1*60*60; tau_y_LKF = 1*60*60; tau_z_LKF = 1*60*60;
            B_LKF = diag([tau_x_LKF^-1, tau_y_LKF^-1, tau_z_LKF^-1]);
            sig_u_LKF = 2.5e-12;
            Qu_LKF = diag([sig_u_LKF^2, sig_u_LKF^2, sig_u_LKF^2]);
            P0_i_LKF = blkdiag(P0_i, Qu_LKF);
            X0_i_LKF = [X0_i; zeros(3,1)];
            x0_i_LKF = zeros(size(X0_i_LKF));
                        % EKF
            tau_x_EKF = 0.825*60*60; tau_y_EKF = 0.825*60*60; tau_z_EKF = 0.825*60*60;
            B_EKF = diag([tau_x_EKF^-1, tau_y_EKF^-1, tau_z_EKF^-1]);
            sig_u_EKF = 2.5e-13;
            Qu_EKF = diag([sig_u_EKF^2, sig_u_EKF^2, sig_u_EKF^2]);
            P0_i_EKF = blkdiag(P0_i, Qu_EKF);

            numLKF = 100;
            LKFiter = 10;
            if numMeas < numLKF
                    % Only run LKF
                LKFRun = runLKF_DMC(X0_i_LKF, x0_i_LKF, P0_i_LKF, B_LKF, Qu_LKF, pConst, scConst, stations, X_days, t_days, tMeas(kStart), numMeas, LKFiter, plotBool_LKF);
                LKFRuns = [LKFRuns; LKFRun];

                    % Accumulate data
                t_Combined = [t_Combined; LKFRun.t_LKF];
                X_Combined = [X_Combined, LKFRun.X_LKF];
                PEst_Combined = [PEst_Combined, LKFRun.LKFOut.PEst];
                prefits_Combined = [prefits_Combined, LKFRun.LKFOut.prefit_res];
                postfits_Combined = [postfits_Combined, LKFRun.LKFOut.postfit_res];
                iterations_Combined = [iterations_Combined, zeros(1,length(LKFRun.t_LKF))];
            else
                    % Run LKF for a few measurements
                LKFRun = runLKF_DMC(X0_i_LKF, x0_i_LKF, P0_i_LKF, B_LKF, Qu_LKF, pConst, scConst, stations, X_days, t_days, tMeas(kStart), numLKF, LKFiter, plotBool_LKF);
                LKFRuns = [LKFRuns; LKFRun];

                    % Accumulate data
                t_Combined = [t_Combined; LKFRun.t_LKF];
                X_Combined = [X_Combined, LKFRun.X_LKF];
                PEst_Combined = [PEst_Combined, LKFRun.LKFOut.PEst];
                prefits_Combined = [prefits_Combined, LKFRun.LKFOut.prefit_res];
                postfits_Combined = [postfits_Combined, LKFRun.LKFOut.postfit_res];
                iterations_Combined = [iterations_Combined, zeros(1,length(LKFRun.t_LKF))];

                    % Run IEKF for the rest
                X0_EKF = X_Combined(:,end); P0_EKF = PEst_Combined{end}; t_EKF = t_Combined(end);
                X0_EKF(8:10) = zeros(3,1); P0_EKF(8:10,8:10) = Qu_EKF;
                EKFRun = runIEKF_DMC(X0_EKF, P0_EKF, B_EKF, Qu_EKF, pConst, scConst, stations, X_days, t_days, t_EKF, tMeas(kEnd_truth), plotBool_EKF);
                EKFRuns = [EKFRuns; EKFRun];

                    % Accumulate data
                t_Combined = [t_Combined; EKFRun.t_EKF];
                X_Combined = [X_Combined, EKFRun.X_EKF];
                PEst_Combined = [PEst_Combined, EKFRun.EKFOut.PEst];
                prefits_Combined = [prefits_Combined, EKFRun.EKFOut.prefit_res];
                postfits_Combined = [postfits_Combined, EKFRun.EKFOut.postfit_res];
                iterations_Combined = [iterations_Combined, EKFRun.EKFOut.iterations];
            end

            % numLKF = 100;
            % LKFiter = 10;
            % if numMeas < numLKF
            %         % Only run LKF
            %     LKFRun = runLKF(X0_i, x0_i, P0_i, pConst, scConst, stations, X_days, t_days, tMeas(kStart), numMeas, LKFiter, plotBool_LKF);
            %     LKFRuns = [LKFRuns; LKFRun];
            % 
            %         % Accumulate data
            %     t_Combined = [t_Combined; LKFRun.t_LKF];
            %     X_Combined = [X_Combined, LKFRun.X_LKF];
            %     PEst_Combined = [PEst_Combined, LKFRun.LKFOut.PEst];
            %     prefits_Combined = [prefits_Combined, LKFRun.LKFOut.prefit_res];
            %     postfits_Combined = [postfits_Combined, LKFRun.LKFOut.postfit_res];
            % else
            %         % Run LKF for a few measurements
            %     LKFRun = runLKF(X0_i, x0_i, P0_i, pConst, scConst, stations, X_days, t_days, tMeas(kStart), numLKF, LKFiter, plotBool_LKF);
            %     LKFRuns = [LKFRuns; LKFRun];
            % 
            %         % Accumulate data
            %     t_Combined = [t_Combined; LKFRun.t_LKF];
            %     X_Combined = [X_Combined, LKFRun.X_LKF];
            %     PEst_Combined = [PEst_Combined, LKFRun.LKFOut.PEst];
            %     prefits_Combined = [prefits_Combined, LKFRun.LKFOut.prefit_res];
            %     postfits_Combined = [postfits_Combined, LKFRun.LKFOut.postfit_res];
            % 
            %         % Run IEKF for the rest
            %     X0_EKF = X_Combined(:,end); P0_EKF = PEst_Combined{end}; t_EKF = t_Combined(end);
            %     EKFRun = runIEKF(X0_EKF, P0_EKF, pConst, scConst, stations, X_days, t_days, t_EKF, tMeas(kEnd_truth), plotBool_EKF);
            %     EKFRuns = [EKFRuns; EKFRun];
            % 
            %         % Accumulate data
            %     t_Combined = [t_Combined; EKFRun.t_EKF];
            %     X_Combined = [X_Combined, EKFRun.X_EKF];
            %     PEst_Combined = [PEst_Combined, EKFRun.EKFOut.PEst];
            %     prefits_Combined = [prefits_Combined, EKFRun.EKFOut.prefit_res];
            %     postfits_Combined = [postfits_Combined, EKFRun.EKFOut.postfit_res];
            % end

            % tau_x = 0.5*60*60; tau_y = 0.5*60*60; tau_z = 0.5*60*60;
            % B = diag([tau_x^-1, tau_y^-1, tau_z^-1]);
            % 
            % sig_u = 7.5e-13;
            % Qu = diag([sig_u^2, sig_u^2, sig_u^2]);
            % 
            % P0_i = blkdiag(P0_i, Qu);
            % 
            % LKFRun = runLKF_DMC(X0_i, x0_i, P0_i, B, Qu, pConst, scConst, stations, X_days(kStart:kEnd_truth,:), t_days(kStart:kEnd_truth), tMeas(kStart), numMeas, 10, plotBool);
            % LKFRuns = [LKFRuns; LKFRun];
            % 
            %     % Accumulate data
            % t_Combined = [t_Combined; LKFRun.t_LKF];
            % X_Combined = [X_Combined, LKFRun.X_LKF];
            % PEst_Combined = [PEst_Combined, LKFRun.LKFOut.PEst];
            % prefits_Combined = [prefits_Combined, LKFRun.LKFOut.prefit_res];
            % postfits_Combined = [postfits_Combined, LKFRun.LKFOut.postfit_res];

            sigma_Combined = [];
            for kk = 1:length(PEst_Combined)
                P = PEst_Combined{kk};
            
                sigPart = [];
                for elem = 1:size(P,1)
                    sigPart = [sigPart, sqrt(P(elem,elem))];
                end
            
                sigma_Combined = [sigma_Combined; sigPart];
            end
            
            % titleText = sprintf("Pre-Fit Residuals over %.0f days", days(k)); 
            % xLabel = "Time [days]"; 
            % yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
            % colors = ['b', 'r'];
            % plotResiduals(t_Combined/(24*60*60), prefits_Combined, titleText, xLabel, yLabel, colors);
            % 
            % titleText = sprintf("Post-Fit Residuals over %.0f days", days(k)); 
            % xLabel = "Time [days]"; 
            % yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
            % colors = ['b', 'r'];
            % plotResiduals(t_Combined/(24*60*60), postfits_Combined, titleText, xLabel, yLabel, colors);
            % 
            % titleText = sprintf("Estimated DMC Accelerations over %.0f days", days(k));
            % xLabel = "Time [days]";
            % yLabel = ["w_X [km/s^2]", "w_Y [km/s^2]", "w_Z [km/s^2]"];
            % colors = ['b', 'b', 'b'];
            % plotW(t_Combined/(24*60*60), X_Combined, titleText, xLabel, yLabel, colors);
            % 
            % %% Part 3: Bplane Implementation
            % fprintf("\n---- Calculating B Plane ----\n")
            % 
            %     % Set up ODE events for sphere of influence
            % opt.Events = @(t,X)SOICheck(t,X,pConst);
            % 
            %     % Define function handle for integration
            % DynFunc = @(t,XPhi)STMEOM_MuSunSRP(t,XPhi,pConst,scConst);
            % 
            %     % Integrate to 3 SOI
            % XPhi_0 = [X_Combined(1:7,end); reshape(eye(7), 49, 1)];
            % P0_Bplane = PEst_Combined{end};
            % P0_Bplane = P0_Bplane(1:7,1:7);
            % tspan_3SOI = (t_Combined(end):1000:300*24*60*60)'; % Add time up to 300 days in sec
            % [t_3SOI, XPhi_3SOI] = ode45(@(t,XPhi)DynFunc(t,XPhi), tspan_3SOI, XPhi_0, opt);
            % Phi_3SOI = reshape(XPhi_3SOI(end, 8:end),7,7);
            % P_3SOI = Phi_3SOI*P0_Bplane*Phi_3SOI';
            % 
            %     % Calculate Bplane without SOI checker
            % [BdotR, BdotT, X_crossing, P_BPlane, STR2ECI, XPhi_BPlane, t_BPlane] = calcBPlane(XPhi_3SOI(end,:)', t_3SOI(end), P_3SOI, pConst, DynFunc, odeset('RelTol',1e-13,'AbsTol',1e-13));
            % 
            %     % Plot BdotR and BdotT location + uncertainty ellipse
            % if k == 1
            %     newFig = true;
            % else
            %     newFig = false;
            % end
            % 
            % boundLevel = 3;
            % titleText = sprintf("Part 3: B Plane Target Estimate from Estimated State");
            % xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; 
            % ellipseLabel = sprintf("%.0f\\sigma B plane target uncertainty, %.3f days of data", boundLevel, t_Combined(end)/(24*60*60));
            % ellipseColor = dayColors(k);
            % BVecLabel = sprintf("\nB Plane Target after %.3f days of data: \nBdotR = %.4e km,\nBdotT = %.4e km", t_Combined(end)/(24*60*60), BdotR, BdotT);
            % plotBPlane(BdotR, BdotT, X_crossing, P_BPlane, STR2ECI, pConst, boundLevel, titleText, xLabel, yLabel, zLabel, ellipseLabel, ellipseColor, BVecLabel, 69, newFig);
            % 
            % figure;
            % hold on; grid on; axis equal;
            % titleText = sprintf("Part 3: 3D Trajectory over %.3f days of data", t_Combined(end)/(24*60*60));
            % title(titleText);    
            % plot3(X_ref(:,1), X_ref(:,2), X_ref(:,3), 'k--')
            % plot3(X_Combined(1,:), X_Combined(2,:), X_Combined(3,:), 'm--')
            % plot3(0,0,0,'b.','MarkerSize', 50);
            % xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
            % view([30 35]);
            % legend("Reference Trajectory", "Estimated Trajectory", "Earth", 'Location', 'best')

        end

            % Ensure 0 iterations show up on plot
        iterations_Combined = iterations_Combined + 0.1;

            % Do plotting and Bplane calculation at specified days
        days = linspace(50,250,5);
        for k = 1:length(days)
                % Find proper indices
            % kStart = find(t_Combined/(24*60*60) >= days(k)-50, 1, 'first');
            kEnd = find(t_Combined/(24*60*60) <= days(k), 1, 'last');

            idx = 1:kEnd;

            if k == length(days)
                    % Plot residuals, W's, and iterations
                titleText = sprintf("Pre-Fit Residuals (Before Iteration) over %.0f days", days(k)); 
                xLabel = "Time [days]"; 
                yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
                colors = ['b', 'r'];
                plotResiduals(t_Combined(idx)/(24*60*60), prefits_Combined(:,idx), titleText, xLabel, yLabel, colors);
            
                titleText = sprintf("Post-Fit Residuals (After Iteration) over %.0f days", days(k)); 
                xLabel = "Time [days]"; 
                yLabel = ["Range Residuals [km]", "Range-Rate Residuals [km/s]"];
                colors = ['b', 'r'];
                plotResiduals(t_Combined(idx)/(24*60*60), postfits_Combined(:,idx), titleText, xLabel, yLabel, colors);
            
                titleText = sprintf("Estimated DMC Accelerations over %.0f days", days(k));
                xLabel = "Time [days]";
                yLabel = ["w_X [km/s^2]", "w_Y [km/s^2]", "w_Z [km/s^2]"];
                colors = ['b', 'b', 'b'];
                plotW(t_Combined(idx)/(24*60*60), X_Combined(:,idx), titleText, xLabel, yLabel, colors);

                titleText = sprintf("IEKF Iterations over %.0f days", days(k));
                xLabel = "Time [days]";
                yLabel = "Measurement Update Iterations";
                colors = 'b';
                figure; 
                sgtitle(titleText)
                subplot(1,2,1)
                    semilogy(t_Combined(idx)/(24*60*60), iterations_Combined, '.', 'Color', colors)
                    hold on; grid on;
                    xlabel(xLabel); ylabel(yLabel);
                subplot(1,2,2)
                    hold on; grid on;
                    histogram(iterations_Combined, 'FaceColor',colors,'Orientation','horizontal')
            end

            %% Part 3: Bplane Implementation
            fprintf("\n---- Calculating B Plane for t = %.0f days ----\n", days(k))
            
                % Set up ODE events for sphere of influence
            opt.Events = @(t,X)SOICheck(t,X,pConst);
            
                % Define function handle for integration
            DynFunc = @(t,XPhi)STMEOM_MuSunSRP(t,XPhi,pConst,scConst);
            
                % Integrate to 3 SOI
            XPhi_0 = [X_Combined(1:7,kEnd); reshape(eye(7), 49, 1)];
            P0_Bplane = PEst_Combined{kEnd};
            P0_Bplane = P0_Bplane(1:7,1:7);
            tspan_3SOI = (t_Combined(kEnd):1000:300*24*60*60)'; % Add time up to 300 days in sec
            [t_3SOI, XPhi_3SOI] = ode45(@(t,XPhi)DynFunc(t,XPhi), tspan_3SOI, XPhi_0, opt);
            Phi_3SOI = reshape(XPhi_3SOI(end, 8:end),7,7);
            P_3SOI = Phi_3SOI*P0_Bplane*Phi_3SOI';
            
                % Calculate Bplane without SOI checker
            [BdotR, BdotT, sig_R, sig_T, sig_RT, X_crossing, P_BPlane, STR2ECI, XPhi_BPlane, t_BPlane] = calcBPlane(XPhi_3SOI(end,:)', t_3SOI(end), P_3SOI, pConst, DynFunc, odeset('RelTol',1e-13,'AbsTol',1e-13));
            
            fprintf("\n\tSigma_R: %.3f, Sigma_T: %.3f, Sigma_RT: %.3e\n", sig_R, sig_T, sig_RT);

                % Plot BdotR and BdotT location + uncertainty ellipse
            if k == 1
                newFig = true;
            else
                newFig = false;
            end
        
            boundLevel = 3;
            titleText = sprintf("Part 3: B Plane Target Estimate from Estimated State");
            xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; 
            ellipseLabel = sprintf("%.0f\\sigma B plane target uncertainty, %.3f days of data", boundLevel, t_Combined(kEnd)/(24*60*60));
            ellipseColor = dayColors(k);
            BVecLabel = sprintf("\nB Plane Target after %.3f days of data: \nBdotR = %.4e km,\nBdotT = %.4e km", t_Combined(kEnd)/(24*60*60), BdotR, BdotT);
            plotBPlane(BdotR, BdotT, X_crossing, P_BPlane, STR2ECI, pConst, boundLevel, titleText, xLabel, yLabel, zLabel, ellipseLabel, ellipseColor, BVecLabel, 69, newFig, part);
        
            % figure;
            % hold on; grid on; axis equal;
            % titleText = sprintf("Part 3: 3D Trajectory over %.3f days of data", t_Combined(end)/(24*60*60));
            % title(titleText);    
            % plot3(X_ref(idx,1), X_ref(idx,2), X_ref(idx,3), 'k--')
            % plot3(X_Combined(1,idx), X_Combined(2,idx), X_Combined(3,idx), 'm--')
            % plot3(0,0,0,'b.','MarkerSize', 50);
            % xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
            % view([30 35]);
            % legend("Reference Trajectory", "Estimated Trajectory", "Earth", 'Location', 'best')

        end

            % Report final estimate of BdotR and BdotT to text file
        result = sprintf("%.3f %.3f", BdotR, BdotT);
        writematrix(result,"Faber_Ian_bplane_overwrite.txt")

    case 4
        fprintf("\nHave a great day!\n")
    otherwise
        fprintf("\nUnknown Option! Have a good day!\n")
end


