function rms = calcResidualRMS(residuals, stations, vis)
% Function that calculates the RMS error for residuals
%   Inputs:
%       - residuals: Filter residials as follows:
%                    prefit residuals (Batch, LKF, and EKF): 
%                       - y_i
%                    postfit residuals (Batch and LKF):
%                       - Batch: y_i - H_i*x_0
%                       - LKF: y_i - Htilde_i*x_i
%       - stations: Stations structure as defined by makeStations.m
%       - vis: Station visibility cell array as defined by
%              processStations.m
%   Outputs:
%       - rms: post-fit residual RMS error
%
%   By: Ian Faber, 02/01/2025
%

rms_part = [];
for kk = 1:size(residuals,2)
    stat = vis{kk}; % Find which station was visible and made the measurement
    rms_part = [rms_part; residuals(:,kk)'*(stations(stat(1)).R{stat(1)}^-1)*residuals(:,kk)];
end
rms = sqrt(sum(rms_part)/(size(residuals,1)*size(residuals,2)));

end