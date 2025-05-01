function [t, Y, R, Xs, vis] = processStations(stations, t_start)
% Function that processes the measurements in the stations struct from
% makeStations.m for multiple station observations at the same time into
% usable arrays
%   Inputs:
%       - stations: stations struct as defined by makeStations.m
%       - t_start: Time after which to start processing measurements
%                  (optional). Defaults to processing all measurements.
%   Outputs:
%       - t: Time vector for measurements
%       - Y: Cell array of measurement vectors, one vector for every time
%            in t:
%            {Y_1, Y_2, Y_3, ..., Y_n}, where n is the number of timesteps
%            where a measurement was made
%       - R: Cell array of measurement covariance matrices, one matrix for
%            every time in t:
%            {R_1, R_2, R_3, ..., R_n}, where n is the number of timesteps
%            where a measurement was made
%       - Xs: Cell array of station states, one matrix for every time in t:
%            {Xs_1, Xs_2, Xs_3, ..., Xs_n}, where n is the number of
%            timesteps where a measurement was made and
%            Xs_i = [Xs_1; Xs_2]
%       - vis: Cell array of station IDs that can see the spacecraft at the
%              time of measurement, one vector for every time in t:
%              {vis_1, vis_2, vis_3, ..., vis_n}, where n is the number of
%              timesteps where a measurement was made
%
%   By: Ian Faber, 01/30/2025
%
    % Initialize outputs
t = [];
Y = {};
R = {};
Xs = {};
vis = {};

    % Pull out all station data, minus elevation angle, and store in
    % trackers for sorting later
measurements = [];
Rs = [];
Xss = [];
viss = [];
for k = 1:length(stations)
    measurements = [measurements; [stations(k).rho, stations(k).rhoDot, stations(k).tMeas]];
    Rs = [Rs; stations(k).R];
    idx = [];
    for kk = 1:length(stations(k).tMeas)
        idx = [idx, find(stations(k).Xs(:,end) == stations(k).tMeas(kk))];
    end
    Xss = [Xss; stations(k).Xs(idx,:)];
    viss = [viss; [stations(k).idx*ones(size(stations(k).tMeas)), stations(k).tMeas]];
end

    % Append time to R matrix tracker
a = measurements(:,3);
Rs = [Rs, mat2cell(a, ones(size(a)))];

    % Sort trackers by time
measurements = sortrows(measurements, 3);
Rs = sortrows(Rs, 2);
Xss = sortrows(Xss,size(Xss,2));
viss = sortrows(viss, 2);

    % Find duplicates - this is when we have multiple measurements at the
    % same time!
time = measurements(:,3);

if exist("t_start", 'var') % If starting at a delayed time, only look at measurements etc. after the start time
    measurements = measurements(time >= t_start,:);
    Rs = Rs(time >= t_start);
    Xss = Xss(time >= t_start,:);
    viss = viss(time >= t_start);
    time = time(time >= t_start);
end

[t, dupeIdx] = unique(time,'stable');

    % Create final observation arrays
for k = 1:length(dupeIdx)
    if k+1 > length(dupeIdx)
        idx = dupeIdx(k);
    else
        idx = dupeIdx(k):dupeIdx(k+1)-1; % Choose indices to create hybrid arrays for multiple measurements
    end

    meas = []; r = []; xs = []; v = [];
    for kk = 1:length(idx)
        meas = [meas; measurements(idx(kk),1:2)'];
        r = blkdiag(r,Rs{idx(kk),1});
        xs = [xs; Xss(idx(kk),1:(size(Xss,2)-1))];
        v = [v; viss(idx(kk),1)];
    end
    Y{k} = meas;
    R{k} = r;
    Xs{k} = xs;
    vis{k} = v;
end

end