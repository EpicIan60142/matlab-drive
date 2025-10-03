function expRange = computeExpectedRange(GPSData, ephData, userPos, pConst)
% Function that calculates expected range at a given time based on
% ephemeris data, estimated user position, and planetary constants
%   Inputs:
%       - GPSData: RINEX GPS data for a given satellite PRN
%       - ephData: Ephemeris data from a satellite broadcast file for a
%                  given PRN, same PRN as the RINEX data
%       - userPos: ECEF user position in meters, organized as [X; Y; Z]
%       - pConst: Planetary constants structure containing the following 
%                 fields:
%                   - wE: rotation rate of earth in rad/s
%                   - c: Speed of light in m/s
%   Outputs:
%       - expRange: Vector of expected range measurements for each time
%                   instant from the RINEX data
%
%   By: Ian Faber, 10/03/2025
%

    % Extract time and calculate WN and TOW
t = GPSData.Time;
gpsEpoch = datetime(1980, 1, 6, 0, 0, 0);
tDiff = t - gpsEpoch;
WN = floor(seconds(tDiff)/(7*24*60*60));
TOW = mod(seconds(tDiff), 7*24*60*60);

    % Extract PRN
PRNRnx = GPSData.SatelliteID(1);
PRNEph = ephData(1,1);

if PRNRnx ~= PRNEph
    fprintf("PRN Mismatch in provided data! Make sure PRNs match...\n")
    expRange = [];
    return;
end

PRN = PRNRnx;

    % Convert ephemeris data to satellite position
[~, ephPos, ~, ~, ~, ~] = eph2pvt2025(ephData, [WN, TOW], PRN);

ephPos = ephPos'; % Make positions column vectors

    % Loop through all measurement times
expRange = [];
for k = 1:length(t)
        % Step 1
    wn_r = WN(k);
    t_r = TOW(k);
    pos_r = ephPos(:, k);

        % Step 2
    R_geo = norm(pos_r - userPos);

        % Step 3
    R_old = R_geo;
    t_t = t_r - R_old/pConst.c;

        % Step 4
    [~, pos_t, ~, ~, ~, ~] = eph2pvt2025(ephData, [wn_r, t_t], PRN);

        % Step 5
    phi = pConst.wE*(t_r - t_t);
    pos_r = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1]*pos_t';

        % Step 6
    R_new = norm(pos_r - userPos);

    while(abs((R_new/R_old)-1) > 1e-2) % Wait for agreement within 1%
        R_old = R_new;

            % Step 3
        t_t = t_r - R_old/pConst.c;

            % Step 4
        [~, pos_t, ~, ~, ~, ~] = eph2pvt2025(ephData, [wn_r, t_t], PRN);
    
            % Step 5
        phi = pConst.wE*(t_r - t_t);
        pos_r = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1]*pos_t';
    
            % Step 6
        R_new = norm(pos_r - userPos);
    end

    expRange = [expRange; R_new];
end


end


