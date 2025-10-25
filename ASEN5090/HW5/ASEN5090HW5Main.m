%% ASEN 5090 HW 5 Main Script
% By: Ian Faber, 10/24/2025

%% Housekeeping
clc; clear; close all;

addpath("../HW3/HW3_Data/");
addpath("../HW3/HW3_Code/");
addpath("../utilities");

%% Setup
    % Truth NIST position
r_NIST_ECEF = [-1288398.567; -4721696.932; 4078625.350]; % m
r_NIST_LLA = ecef2lla(r_NIST_ECEF');

    % L1 and L2 frequencies
f1 = 1575.42e6; % Hz
f2 = 1227.60e6; % Hz

    % Planetary constants structure
pConst.c = 299792458; % m/s
pConst.wE = 7.2921151467e-5; % rad/s

    % tropospheric zenith delay and elevation cutoff
zd = 2; % m
elevationCutoff = 5; % deg    

    % Read RINEX data
rinexData = rinexread("NIST00USA_R_20252230000_01D_30S_MO.rnx");
GPSDataFull = rinexData.GPS;
rinexTimes = unique(GPSDataFull.Time);

    % Read ephemeris data
ephData = read_clean_GPSbroadcast("brdc2230.25n", true);

    % Declare PRNs of interest
PRN = [1, 3, 4, 6, 9, 16, 26, 28, 31, 32];

    % Extract only those PRNs
PRNRnxIdx = logical(sum(GPSDataFull.SatelliteID == PRN, 2));
GPSData = GPSDataFull(PRNRnxIdx, :);

PRNEphIdx = logical(sum(ephData(:,1) == PRN, 2));
ephData = ephData(PRNEphIdx, :);

%% Problem 1: Data gathering
    % Part a. Find distance to initial guess
format long
r_0_LLA = [40, -105, 1600];
r_0_ECEF = lla2ecef(r_0_LLA)'
% r_0_ECEF = lla2ecef(r_NIST_LLA)';

dr_0_ECEF = r_NIST_ECEF - r_0_ECEF
dist_0 = norm(dr_0_ECEF)

    % Part b. Find PIF, R1, El, Az, Bsv, Rel, and Trop for each PRN in the
    % RINEX file
        % Find indices corresponding to the first measurement epoch
idxFirstEpoch = GPSData.Time == rinexTimes(1);

        % Pull out C1C and C2W for the first epoch
C1C_first = GPSData.C1C(idxFirstEpoch);
C2W_first = GPSData.C2W(idxFirstEpoch);

        % Calculate PIF
[PIF_first, ~] = ionocorr(C1C_first, f1, C2W_first, f2)

        % Calculate R1
[~, ~, R1, posT] = computeExpectedRange(GPSData(idxFirstEpoch, :), ephData, r_0_ECEF, pConst);
R1_first = cell2mat(R1)
posT = cell2mat(posT')

        % Calculate Elevation and Azimuth at T_transmit
[satAz_first, satEl_first, ~] = compute_azelrange(r_0_ECEF, posT)

        % Calculate clock corrections for each PRN at the first epoch
clkCorr_first = [];
relCorr_first = [];
for k = 1:length(PRN)
            % Calculate times and convert to Week Number of Time Of Week
    t = GPSData.Time(PRN(k));
    gpsEpoch = datetime(1980, 1, 6, 0, 0, 0);
    tDiff = t - gpsEpoch;
    WN = floor(seconds(tDiff)/(7*24*60*60));
    TOW = mod(seconds(tDiff), 7*24*60*60);

            % Calculate clock and relativity corrections
    try
        [~, ~, ~, clkCorr, relCorr, ~, ~] = eph2pvtWithRel(ephData, [WN, TOW], PRN(k));
    catch
        fprintf("\n\tBad ephemeris data for PRN %.0f! Can't calculate clock corrections\n", PRN(k));
        continue;
    end

            % Append data to vector
    clkCorr_first = [clkCorr_first; clkCorr];
    relCorr_first = [relCorr_first; relCorr];
end
clkCorr_first
relCorr_first

        % Calculate tropospheric model
[tropoCorr_first, ~] = tropomodel(zd, deg2rad(satEl_first), elevationCutoff)

%% Problem 2: First Epoch point solution
    % Part a. Calculate and print the prefit residuals for each satellite
        % PRN 32 doesn't have a valid PIF, look at rows 1:9
rho = PIF_first(1:9) + clkCorr_first + relCorr_first - tropoCorr_first; % Measurement and corrections
rho_0 = R1_first + 0; % Expected range and bias
prefit_first = rho - rho_0;

    % Part b. Calculate and print the measurement sensitivity matrix
G = [-(posT'-r_0_ECEF')./R1_first, ones(size(posT,2), 1)]

    % Part c. Calculate least squares solution for the first epoch, without
    % iterating
dx = ((G'*G)^-1)*G'*prefit_first

    % Part d. Update your estimate of the position and compare to true NIST
    % location
        % Extract updates
dr = dx(1:3); db = dx(4);

        % Apply updates
r_1 = r_0_ECEF + dr
b_1 = 0 + db

        % Find distance to truth
dr_1 = r_NIST_ECEF - r_1
dist_1 = norm(dr_1)

    % Part e. Iterate until norm(dx) < 0.01
iter = 0;
epsilon = 0.01;
while(norm(dx) >= epsilon)
        % Recompute expected range
    [~,~,R1,posT] = computeExpectedRange(GPSData(idxFirstEpoch, :), ephData, r_1, pConst);
    R1 = cell2mat(R1);
    posT = cell2mat(posT');

        % Recompute corrections
            % New elevation for tropospheric delay
    [~, satEl, ~] = compute_azelrange(r_1, posT);
            % Clock corrections
    clkCorr_1 = [];
    relCorr_1 = [];
    for k = 1:length(PRN)
                % Calculate times and convert to Week Number of Time Of Week
        t = GPSData.Time(PRN(k));
        gpsEpoch = datetime(1980, 1, 6, 0, 0, 0);
        tDiff = t - gpsEpoch;
        WN = floor(seconds(tDiff)/(7*24*60*60));
        TOW = mod(seconds(tDiff), 7*24*60*60);
    
                % Calculate clock and relativity corrections
        try
            [~, ~, ~, clkCorr, relCorr, ~, ~] = eph2pvtWithRel(ephData, [WN, TOW], PRN(k));
        catch
            fprintf("\n\tBad ephemeris data for PRN %.0f! Can't calculate clock corrections\n", PRN(k));
            continue;
        end
    
                % Append data to vector
        clkCorr_1 = [clkCorr_1; clkCorr];
        relCorr_1 = [relCorr_1; relCorr];
    end
            % Troposphere
    tropoCorr_1 = tropomodel(zd, satEl, elevationCutoff);

        % Compute prefit residual
    rho = PIF_first(1:9) + clkCorr_1 + relCorr_1 - tropoCorr_1;
    rho_0 = R1 + b_1;
    dRho = rho - rho_0;

        % Compute sensitivity matrix
    G = [-(posT'-r_0_ECEF')./R1_first, ones(size(posT,2), 1)];

        % Calculate correction
    dx = ((G'*G)^-1)*G'*dRho;

        % Apply correction
    r_1 = r_1 + dx(1:3);
    b_1 = b_1 + dx(4);
    iter = iter + 1;

    fprintf("\nIteration %.0f - Correction: %.4f m, Position: [%.2f, %.2f, %.2f]\n", iter, norm(dx), r_1);

end

    % Part f. Compute and display the postfit residuals and the final error
b_f = b_1

r_f = r_1
dr_f = r_NIST_ECEF - r_f

postfits = dRho - G*dx

    % Part g. Report dr_f in ENU coordinates and calculate HDOP and VDOP
        % Get DCM to ENU from ECEF
C_ECEF2ENU = ECEF2ENU(r_NIST_LLA(1), r_NIST_LLA(2));
        % Convert coordinates
dr_f_ENU = C_ECEF2ENU*dr_f

        % Calculate H in ECEF, then convert to ENU
H_ECEF = (G'*G)^-1;
C = blkdiag(C_ECEF2ENU, 1); % Extend the DCM with 1, time doesn't need to change coordinates
H_ENU = C*H_ECEF*C';

HDOP = sqrt(H_ENU(1,1) + H_ENU(2,2))
VDOP = sqrt(H_ENU(3,3))

%% Problem 3: Least squares point solutions for all epochs
    % Use all PRNs
GPSData = GPSDataFull;

    % Preallocate storage vectors
times = [];
positions = [];
clockBiases = [];
positionErrors = [];
HDOPs = [];
VDOPs = [];
nSolSats = [];
postfits = [];
postfitTimes = [];

    % Calculate ECEF to ENU DCM
C_ECEF2ENU = ECEF2ENU(r_NIST_LLA(1), r_NIST_LLA(2));

    % Loop over every time epoch and solve the least squares problem
maxIter = 1000;
for k = 1:length(rinexTimes)
        % Pull out time index
    timeIdx = GPSData.Time == rinexTimes(k);

        % Calculate PIF using C1C and C2W
    C1C = GPSData.C1C(timeIdx);
    C2W = GPSData.C2W(timeIdx);
    [PIF_start, ~] = ionocorr(C1C, f1, C2W, f2);

    if(length(PIF_start) < 4)
        fprintf("\n\tNot enough satellites visible for a solution\n")
        continue;
    end

        % Start with an initial receiver location and bias guess
    rHat = r_0_ECEF;
    bHat = 0;

        % Loop until correction vector is less than 1 cm
    dx = 9999*ones(4,1);
    iter = 0;
    breakEarly = false;
    while(norm(dx) > 0.01 && iter < maxIter)
            % Compute expected range
        [~, ~, R1, posT] = computeExpectedRange(GPSData(timeIdx, :), ephData, rHat, pConst);
        R1 = cell2mat(R1);
        posT = cell2mat(posT');
    
            % Compute clock corrections
        PRN = GPSData(timeIdx,:).SatelliteID;
        PRNInvalid = false(size(PRN));
        clkCorr = [];
        relCorr = [];
        for kk = 1:length(PRN)
                    % Calculate times and convert to Week Number of Time Of Week
            t = GPSData.Time(PRN(kk));
            gpsEpoch = datetime(1980, 1, 6, 0, 0, 0);
            tDiff = t - gpsEpoch;
            WN = floor(seconds(tDiff)/(7*24*60*60));
            TOW = mod(seconds(tDiff), 7*24*60*60);
        
                    % Calculate clock and relativity corrections
            try
                [~, ~, ~, clkCorrPRN, relCorrPRN, ~, ~] = eph2pvtWithRel(ephData, [WN, TOW], PRN(kk));
            catch
                PRNInvalid(kk) = true;
                continue;
            end
        
                    % Append data to vector
            clkCorr = [clkCorr; clkCorrPRN];
            relCorr = [relCorr; relCorrPRN];
        end
    
            % Update satellite position and expected range validity
        % R1 = R1(~PRNInvalid);
        % posT = posT(:,~PRNInvalid);

            % Troposphere
        [~, satEl, ~] = compute_azelrange(rHat, posT);
        [tropoCorr, idxValid] = tropomodel(zd, satEl, elevationCutoff);

            % Update vector validity
        PIF_goodEph = PIF_start(~PRNInvalid); % Remove bad ephemeris data
        if length(PIF_goodEph) < 4
            fprintf("\n\tNot enough satellites visible for solution\n")
            breakEarly = true;
            break;
        end
        PIF = PIF_goodEph(idxValid); % Remove low elevations
        clkCorr = clkCorr(idxValid);
        relCorr = relCorr(idxValid);
        tropoCorr = tropoCorr(idxValid);
        R1 = R1(idxValid);
        posT = posT(:, idxValid);
        
            % Compute prefits
        rho = PIF + clkCorr + relCorr - tropoCorr;
        rho_0 = R1 + bHat;
        dRho = rho - rho_0;

            % Compute measurement sensitivity matrix
        G = [-(posT'-rHat')./R1, ones(size(posT,2), 1)];

            % Solve for correction
        dx = ((G'*G)^-1)*G'*dRho;

            % Apply correction
        rHat = rHat + dx(1:3);
        bHat = bHat + dx(4);
        iter = iter + 1;
    end

    if breakEarly
        continue;
    end

        % Calculate final position error in ENU coordinates
    rError_ECEF = r_NIST_ECEF - rHat;
    rError_ENU = C_ECEF2ENU*rError_ECEF;

        % Calculate HDOP and VDOP
    H_ECEF = (G'*G)^-1;
    C = blkdiag(C_ECEF2ENU, 1); % Extend the DCM with 1, time doesn't need to change coordinates
    H_ENU = C*H_ECEF*C';
    
    HDOP = sqrt(H_ENU(1,1) + H_ENU(2,2));
    VDOP = sqrt(H_ENU(3,3));
    if ~isreal(HDOP) || ~isreal(VDOP)
        continue;
    end

        % Pull out number of satellites used for solution
    nSats = length(dRho);

        % Calculate postfits
    postfit = dRho - G*dx;
    PRN = PRN(~PRNInvalid);
    PRN = PRN(idxValid);
    postTime = repmat(rinexTimes(k), length(postfit), 1);

        % Accumulate data
    times = [times, rinexTimes(k)];
    positions = [positions, rHat];
    clockBiases = [clockBiases, bHat];
    positionErrors = [positionErrors, rError_ENU];
    HDOPs = [HDOPs, HDOP];
    VDOPs = [VDOPs, VDOP];
    nSolSats = [nSolSats, nSats];
    postfits = [postfits; [postfit, PRN]];
    postfitTimes = [postfitTimes; postTime];
end

    % Plot position errors and clock bias
figure; tl = tiledlayout(4,1); ax = [];
title(tl, "ENU Position Errors and Clock Bias")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(times, positionErrors(1,:), 'b.');
    ylabel("\Delta x_E [m]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(times, positionErrors(2,:), 'b.');
    ylabel("\Delta x_N [m]");    
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(times, positionErrors(3,:), 'b.');
    ylabel("\Delta x_U [m]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(times, clockBiases, 'r.');
    xlabel("Time [GPS]"); ylabel("Estimated clock Bias [m]");
linkaxes(ax,'x');

    % Plot number of satellites used and HDOP/VDOP
figure; tl = tiledlayout(3,1); ax = [];
title(tl, "Number of Satellites Used and HDOP/VDOP")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(times, nSolSats, 'b.');
    ylabel("Number of Satellites used");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(times, HDOPs, 'r.');
    ylabel("HDOP [m]");    
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(times, VDOPs, 'k.');
    xlabel("Time [GPS]"); ylabel("VDOP [m]");

    % Plot postfit residuals
figure;
hold on; grid on;
title("Postfits by PRN")
scatter(postfitTimes, postfits(:,1), 5, postfits(:,2), 'filled')
xlabel("Time [GPS]"); ylabel("Postfit Residual [m]");
cBar = colorbar; cBar.Label.String = "PRN";
colormap("cool")




