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
scConst = getSCConst();

    % Initialize stations struct
stations = makeStations(pConst, 2);

    % ode45 options
opt = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

%% Part 1: Dynamics Implementation
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
title("Part 1: Truth orbit comparison until 200 days - Sun origin frame")
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
title("Part 1: Truth orbit comparison until 200 days - Earth origin frame")
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

    % Propagate station states and generate clean measurements
for k = 1:length(tMeas)
    dTheta = pConst.wEarth*t_test(k);
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
[BdotR_truth, BdotT_truth, X_crossing_truth, P_BPlane, STR2ECI, XPhi_BPlane, t_BPlane] = calcBPlane(XPhi_3SOI(end,:)', t_3SOI(end), P_3SOI, pConst, DynFunc, odeset('RelTol',1e-13,'AbsTol',1e-13));

    % Plot BdotR and BdotT location + uncertainty ellipse
boundLevel = 3;
titleText = sprintf("Part 1: B Plane Target Estimate from Truth Data");
xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; 
ellipseLabel = sprintf("\\pm %.0f\\sigma B plane target uncertainty, dummy P_0", boundLevel);
ellipseColor = 'b';
BVecLabel = sprintf("B Plane Target after propagating truth data: \nBdotR = %.4e km,\nBdotT = %.4e km", BdotR_truth, BdotT_truth);
plotBPlane(BdotR_truth, BdotT_truth, X_crossing_truth, P_BPlane, STR2ECI, pConst, boundLevel, titleText, xLabel, yLabel, zLabel, ellipseLabel, ellipseColor, BVecLabel, 420, true);

figure(420)
nexttile(1)
    hold on;
    idx = length(t_BPlane); offset = 1000;
    orbit = plot3(XPhi_BPlane(idx-offset:end,1), XPhi_BPlane(idx-offset:end,2), XPhi_BPlane(idx-offset:end,3), 'm--', 'DisplayName', "Modeled Orbit");

% return

%% Part 2: Estimate State with Known Target and Models with varying amounts of data
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
for k = 1:length(days)
    fprintf("\n--- Analyzing %.0f days of data ---\n", days(k));

        % Create new truth data for the specified number of days
    if k == 1
        kStart = 1;   
    else
        kStart = kEnd + 1;
    end
    kEnd_truth = find(tMeas <= days(k)*24*60*60, 1, 'last'); 

    XPhi_0 = [truthTraj.Xt_50(1,1:7)'; reshape(eye(length(X0)),length(X0)^2,1)];
    [t_days, XPhi_days] = ode45(@(t,X)STMEOM_MuSunSRP(t,X,pConst,scConst), tMeas(1:kEnd_truth), XPhi_0, opt);
    X_days = XPhi_days(:,1:7);

        % Properly initialize filters
    if k == 1
        X0_i = X0;
        x0_i = x0;
        P0_i = P0;
    else
        X0_i = X_days(kStart, :)';
        x0_i = x0;
        P0_i = 0.5*P0; % Start new segment with high uncertainty to avoid saturation
    end

        % Determine cutoffs between LKF and EKF, i.e. when a large measurement
        % gap happens
    dt = diff(t_days); % Find time between measurements
    med_dt = median(dt); % Find median time gap
    kSec = [find(dt > med_dt); kEnd_truth]; % Pull out when time gaps are larger than the median
    
        % Configure plotting
            % Only plot residuals and state errors w/ bounds
    plotBool = [false; false; false; false; false; false; false; false];     
    
        % Run LKF/EKF, initializing with LKF after a large gap
    % LKFRuns = [];
    % EKFRuns = [];
    % t_Combined = [];
    % X_Combined = [];
    % PEst_Combined = {};
    % prefits_Combined = [];
    % postfits_Combined = [];
    % X0_i = X0;
    % x0_i = x0;
    % P0_i = P0;
    n = 100;
    % kStart = 1;
    if k == 1
        start = 1;
    else
        start = sec + 1;
    end

    for sec = start:length(kSec)
        kEnd = kSec(sec);
    
        if tMeas(kStart) == tMeas(kEnd) % Only 1 measurement
            continue;
        end

        fprintf("\n---- Filtering from t = %.3f to t  = %.3f ----\n", tMeas(kStart), tMeas(kEnd));
    
            % Determine when to turn on EKF
        if kEnd - kStart < n % Use only LKF
            numMeas = kEnd - kStart;
        else % Initialize with n LKF measurements, then switch to EKF
            numMeas = n;
        end
    
        % if k == 3
        %     buh = 1;
        % end
    
            % Run LKF
        LKFRun = runLKF(X0_i, x0_i, P0_i, pConst, scConst, stations, X_days(kStart:kEnd,:), t_days(kStart:kEnd), tMeas(kStart), numMeas, 10, plotBool);
        LKFRuns = [LKFRuns; LKFRun];
        P0_i = LKFRun.LKFOut.PEst{end};
    
            % Accumulate data
        t_Combined = [t_Combined; LKFRun.t_LKF];
        X_Combined = [X_Combined, LKFRun.X_LKF];
        PEst_Combined = [PEst_Combined, LKFRun.LKFOut.PEst];
        prefits_Combined = [prefits_Combined, LKFRun.LKFOut.prefit_res];
        postfits_Combined = [postfits_Combined, LKFRun.LKFOut.postfit_res];
    
            % Run LKF if we have enough measurments to start
        if numMeas == n
            t_start = LKFRun.t_LKF(end); t_end = tMeas(kEnd);

            if t_start == t_end % Only 1 measurement
                continue;
            end
            
            P_start = LKFRun.LKFOut.PEst{end}; X_start = LKFRun.X_LKF(:,end);
    
            EKFRun = runEKF(X_start, P_start, pConst, scConst, stations, X_days(kStart:kEnd,:), t_days(kStart:kEnd), t_start, t_end, plotBool);
            EKFRuns = [EKFRuns; EKFRun];
    
            P0_i = EKFRun.EKFOut.PEst{end};
    
                    % Accumulate data
            t_Combined = [t_Combined; EKFRun.t_EKF];
            X_Combined = [X_Combined, EKFRun.X_EKF];
            PEst_Combined = [PEst_Combined, EKFRun.EKFOut.PEst];
            prefits_Combined = [prefits_Combined, EKFRun.EKFOut.prefit_res];
            postfits_Combined = [postfits_Combined, EKFRun.EKFOut.postfit_res];
        end
    
            % Reset for next set of measurements
        if sec < length(kSec)
                % Update starting index
            kStart = kEnd + 1;
    
                % Get starting covariance for next time
            Phi = eye(length(X0));
            XPhi_0 = [X_Combined(:,end); reshape(Phi, length(X0)^2,1)];
            [~,XPhi] = ode45(@(t,XPhi)STMEOM_MuSunSRP(t,XPhi,pConst,scConst), [t_Combined(end), t_days(kStart)], XPhi_0, opt);
            Phi = reshape(XPhi(end,(length(X0)+1):end), length(X0), length(X0));
    
            P0_i = Phi*PEst_Combined{end}*Phi';
    
                % Update X0 and x0
            X0_i = X_days(kStart,:)';
            x0_i = zeros(size(X0));
        end
    end
    
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
    
    titleText = sprintf("Part 2: Nominal vs. Estimated State over %.3f days", t_Combined(end)/(24*60*60));
    xLabel = "Time [sec]";
    yLabel = ["X [km]", "Y [km]", "Z [km]", ...
              "Xdot [km/s]", "Ydot [km/s]", "Zdot [km/s]", ...
              "C_R [n.d.]"];
    plotStateError(t_Combined, X_Combined', [], [], [], titleText, xLabel, yLabel);
    for kk = 1:size(X_days,2)
        nt = nexttile(kk);
            hold on; grid on;
            plot(t_days, X_days(:,kk),'k--','DisplayName',"Nominal Orbit");
    end
    
    %% Part 2: Bplane Implementation
    fprintf("\n---- Calculating B Plane ----\n")
    
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
    [BdotR, BdotT, X_crossing, P_BPlane, STR2ECI, XPhi_BPlane, t_BPlane] = calcBPlane(XPhi_3SOI(end,:)', t_3SOI(end), P_3SOI, pConst, DynFunc, odeset('RelTol',1e-13,'AbsTol',1e-13));
    
        % Plot BdotR and BdotT location + uncertainty ellipse
    if k == 1
        newFig = true;
    else
        newFig = false;
    end

    boundLevel = 3;
    titleText = sprintf("Part 2: B Plane Target Estimate from Estimated State");
    xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; 
    ellipseLabel = sprintf("\\pm %.0f\\sigma B plane target uncertainty, %.3f days of data", boundLevel, t_Combined(end)/(24*60*60));
    ellipseColor = dayColors(k);
    BVecLabel = sprintf("\nB Plane Target after %.3f days of data: \n\tBdotR = %.4e km,\n\tBdotT = %.4e km", t_Combined(end)/(24*60*60), BdotR, BdotT);
    plotBPlane(BdotR, BdotT, X_crossing, P_BPlane, STR2ECI, pConst, boundLevel, titleText, xLabel, yLabel, zLabel, ellipseLabel, ellipseColor, BVecLabel, 69, newFig);
    
    % figure(69)
    % hold on;
    % idx = length(t_BPlane); offset = 200;
    % orbit = plot3(XPhi_BPlane(idx-offset:end,1), XPhi_BPlane(idx-offset:end,2), XPhi_BPlane(idx-offset:end,3), 'm--', 'DisplayName', "Propagated Orbit");

end

% titleText = sprintf("Part 2: B Plane Target Estimate from Truth Data");
xLabel = "X [km]"; yLabel = "Y [km]"; zLabel = "Z [km]"; 
ellipseLabel = sprintf("\\pm %.0f\\sigma B plane target uncertainty, 200 days of data", boundLevel);
ellipseColor = 'g';
BVecLabel = sprintf("\nB Plane Target from truth data: \n\tBdotR = %.4e km,\n\tBdotT = %.4e km", BdotR_truth, BdotT_truth);
plotBPlane(BdotR_truth, BdotT_truth, X_crossing, P_BPlane, STR2ECI, pConst, boundLevel, titleText, xLabel, yLabel, zLabel, ellipseLabel, ellipseColor, BVecLabel, 69, false);


