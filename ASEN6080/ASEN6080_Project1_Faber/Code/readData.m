function stations = readData(stations,data)
% Parses provided data from a text file into a valid stations struct format
%   - Inputs: 
%       - stations: Stations struct as defined by makeStations.m
%       - data: Data read from a text file
%   - Outputs:
%       - stations: Stations struct with properly formatted measurement
%                   data
%
%   By: Ian Faber, 02/11/2025
%

tMeas = data(:,1);
id = data(:,2);
rho = data(:,3)/1000;
rhoDot = data(:,4)/1000;

ids = [];
for k = 1:length(stations)
    ids = [ids; stations(k).id];
end

for k = 1:length(tMeas)
    statID = id(k);
    idx = find(ids == statID);
    
    stations(idx).tMeas = [stations(idx).tMeas; tMeas(k)];
    stations(idx).rho = [stations(idx).rho; rho(k)];
    stations(idx).rhoDot = [stations(idx).rhoDot; rhoDot(k)];
    stations(idx).R = [stations(idx).R; {diag([stations(idx).sigRho^2, stations(idx).sigRhoDot^2])}];
end


end