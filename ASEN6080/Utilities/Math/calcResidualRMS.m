function rms = calcResidualRMS(residuals, stations, vis, measInclude)
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
%       - measInclude: boolean array indicating which measurements to
%                      include in the residual calculation. To process both
%                      range and range rate, pass in true(1,2) for example.
%   Outputs:
%       - rms: post-fit residual RMS error
%
%   By: Ian Faber, 02/01/2025
%

rms_part = [];
for kk = 1:size(residuals,2)
    try
        stat = vis{kk}; % Find which station was visible and made the measurement
    catch
        stat = vis(kk);
    end
    R = [];
    measCov = stations(stat(1)).R{stat(1)};
    if measInclude(1)
        R = blkdiag(R, measCov(1,1));
    end
    if measInclude(2)
        R = blkdiag(R, measCov(2,2));
    end
    rms_part = [rms_part; residuals(:,kk)'*(R^-1)*residuals(:,kk)];
end
rms = sqrt(sum(rms_part)/(size(residuals,1)*size(residuals,2)));

end