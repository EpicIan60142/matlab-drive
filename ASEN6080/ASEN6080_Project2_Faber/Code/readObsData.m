function [tMeas, stations] = readObsData(stations,data)
% Parses provided data from a text file into a valid stations struct format
% for Project 2
%   - Inputs: 
%       - stations: Stations struct as defined by makeStations.m
%       - data: Data read from a text file
%   - Outputs:
%       - stations: Stations struct with properly formatted measurement
%                   data
%
%   By: Ian Faber, 04/06/2025
%

tMeas = data(:,1);
range = data(:,2:4);
rangeRate = data(:,5:7);

% ids = [];
% for k = 1:length(stations)
%     ids = [ids; stations(k).id];
% end

for k = 1:length(tMeas)
    % statID = id(k);
    % idx = find(ids == statID);

    idx = ~isnan(range(k,:)); % Only assign data if it's not NaN

    stations(idx).tMeas = [stations(idx).tMeas; tMeas(k)];
    stations(idx).rho = [stations(idx).rho; range(k,idx)];
    stations(idx).rhoDot = [stations(idx).rhoDot; rangeRate(k,idx)];
    stations(idx).R = [stations(idx).R; {diag([stations(idx).sigRho^2, stations(idx).sigRhoDot^2])}];
end


end