function [expRange, posT, satRanges, satPosTs] = computeExpectedRange(GPSData, ephData, userPos, pConst)
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
%                   instant from the RINEX data for the last PRN in
%                   GPSData
%       - posT: Vector of positions of the satellite at transmit time for
%               each time instant from the RINEX data for the last PRN in
%               GPSData
%       - satRanges: Cell array of expected ranges for each PRN in GPSData
%       - satPosTs: Cell array of PRN positions at transmit time for each
%                   PRN in GPSData
%
%   By: Ian Faber, 10/03/2025
%
    % Extract PRN
PRNRnx = unique(GPSData.SatelliteID);
PRNEph = unique(ephData(:,1));

if length(PRNRnx) < length(PRNEph)
    PRN = PRNRnx;
else
    PRN = PRNEph;
end

    % Loop over each PRN if more than 1 is provided
satRanges = {};
satPosTs = {};
for kk = 1:length(PRN)
        % Find PRN index
    PRNIdx = GPSData.SatelliteID == PRN(kk);

        % Extract time and calculate WN and TOW
    t = GPSData.Time(PRNIdx);
    gpsEpoch = datetime(1980, 1, 6, 0, 0, 0);
    tDiff = t - gpsEpoch;
    WN = floor(seconds(tDiff)/(7*24*60*60));
    TOW = mod(seconds(tDiff), 7*24*60*60);

        % Convert ephemeris data to satellite position
    try
        [~, ephPos, ~, ~, ~, ~] = eph2pvt2025(ephData, [WN, TOW], PRN(kk));
    catch
        fprintf("\n\tBad ephemeris data for PRN %.0f! Can't calculate expected range\n", PRN(kk));
        continue;
    end

    ephPos = ephPos'; % Make positions column vectors
    
        % Loop through all measurement times
    expRange = [];
    posT = [];
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
        [~, pos_t, ~, ~, ~, ~] = eph2pvt2025(ephData, [wn_r, t_t], PRN(kk));
    
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
            [~, pos_t, ~, ~, ~, ~] = eph2pvt2025(ephData, [wn_r, t_t], PRN(kk));
        
                % Step 5
            phi = pConst.wE*(t_r - t_t);
            pos_r = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1]*pos_t';
        
                % Step 6
            R_new = norm(pos_r - userPos);
        end
    
        expRange = [expRange; R_new];
        posT = [posT, pos_t'];
    end
    satRanges = [satRanges; {expRange}];
    satPosTs = [satPosTs; {posT}];
end



end
