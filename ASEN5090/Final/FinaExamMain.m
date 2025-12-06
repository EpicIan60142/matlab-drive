%% ASEN 5090 Final Exam Main Script
% By: Ian Faber, 12/05/2025

%% Housekeeping
clc; clear; close all;

addpath("../utilities/")

%% Problem 2 Setup
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
rinexData = rinexread("MYST00000_R_20252570000_SINGLE_EPOCH.rnx");
GPSDataFull = rinexData.GPS;
rinexTimes = unique(GPSDataFull.Time);

    % Read ephemeris data
ephDataFull = read_clean_GPSbroadcast("brdc2570.25n", true);

    % Declare PRNs of interest
PRN = [3, 4, 6, 9, 11, 12, 16, 25, 26, 28, 29, 31, 32];

    % Extract only those PRNs
PRNRnxIdx = logical(sum(GPSDataFull.SatelliteID == PRN, 2));
GPSData = GPSDataFull(PRNRnxIdx, :);

PRNEphIdx = logical(sum(ephDataFull(:,1) == PRN, 2));
ephData = ephDataFull(PRNEphIdx, :);

%% Data Gathering
    % Guess an initial position
r_0_LLA = [70, -35, 0];
r_0_ECEF = lla2ecef(r_0_LLA)'

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
    PRNIdx = GPSData.SatelliteID == PRN(k);
    t = GPSData.Time(idxFirstEpoch & PRNIdx);
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

idxNaN = isnan(tropoCorr_first);

%% Find Epoch point solution
    % Part a. Calculate and print the prefit residuals for each satellite
rho = PIF_first(~idxNaN) + clkCorr_first(~idxNaN) + relCorr_first(~idxNaN) - tropoCorr_first(~idxNaN); % Measurement and corrections
rho_0 = R1_first(~idxNaN) + 0; % Expected range and bias
prefit_first = rho - rho_0

    % Part b. Calculate and print the measurement sensitivity matrix
G = [-(posT(:, ~idxNaN)'-r_0_ECEF')./R1_first(~idxNaN), ones(size(posT(:,~idxNaN),2), 1)]

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
        PRNIdx = GPSData.SatelliteID == PRN(k);
        t = GPSData.Time(idxFirstEpoch & PRNIdx);
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
    tropoCorr_1 = tropomodel(zd, deg2rad(satEl), elevationCutoff);
    idxNaN = isnan(tropoCorr_1);

        % Compute prefit residual
    rho = PIF_first(~idxNaN) + clkCorr_1(~idxNaN) + relCorr_1(~idxNaN) - tropoCorr_1(~idxNaN);
    rho_0 = R1(~idxNaN) + b_1;
    dRho = rho - rho_0;

        % Compute sensitivity matrix
    G = [-(posT(:,~idxNaN)'-r_0_ECEF')./R1_first(~idxNaN), ones(size(posT(:,~idxNaN),2), 1)];

        % Calculate correction
    dx = ((G'*G)^-1)*G'*dRho;

        % Apply correction
    r_1 = r_1 + dx(1:3);
    b_1 = b_1 + dx(4);
    iter = iter + 1;

    fprintf("\nIteration %.0f - Correction: %.4f m, Position: [%.2f, %.2f, %.2f], Clock bias: %.2f\n", iter, norm(dx), r_1, b_1);
end

    % Part f. Compute and display the postfit residuals and the final error
b_f = b_1

r_f = r_1

r_LLA = ecef2lla(r_f')

postfits = dRho - G*dx

    % Part g. Report dr_f in ENU coordinates and calculate HDOP and VDOP
        % Get DCM to ENU from ECEF
C_ECEF2ENU = ECEF2ENU(r_0_LLA(1), r_0_LLA(2));

        % Calculate H in ECEF, then convert to ENU
H_ECEF = (G'*G)^-1;
C = blkdiag(C_ECEF2ENU, 1); % Extend the DCM with 1, time doesn't need to change coordinates
H_ENU = C*H_ECEF*C';

HDOP = sqrt(H_ENU(1,1) + H_ENU(2,2))
VDOP = sqrt(H_ENU(3,3))





