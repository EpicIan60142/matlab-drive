function rms = calcPostFitRMS(epsilon, stations, vis)
% Function that calculates the RMS error for postfit residuals
%   Inputs:
%       - epsilon: Filter observation errors 
%                  (y_i - H_i*x_i or y_i - Htilde_i*x_i)
%       - stations: Stations structure as defined by makeStations.m
%       - vis: Station visibility cell array as defined by
%              processStations.m
%   Outputs:
%       - rms: post-fit residual RMS error
%
%   By: Ian Faber, 02/01/2025
%

rms_part = [];
for kk = 1:size(epsilon,2)
    stat = vis{kk}; % Find which station was visible and made the measurement
    rms_part = [rms_part; epsilon(:,kk)'*(stations(stat(1)).R{stat(1)}^-1)*epsilon(:,kk)];
end
rms = sqrt(sum(rms_part)/(size(epsilon,1)*size(epsilon,2)));

end