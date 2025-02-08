function [rms_comp, rms_full] = calcStateErrorRMS(stateError)
% Function that calculates the component-wise and state-wise RMS errors 
% after filtering
%   Inputs:
%       - stateError: State error matrix - one row per timestep and one 
%                     column per state
%   Outputs:
%       - rms_comp: Component-wise RMS error
%       - rms_full: State-wise RMS error
%
%   By: Ian Faber, 02/01/2025
%

    % Component-wise RMS
rms_comp = [];
rms_part = stateError.^2;
for k = 1:size(rms_part,2)
    rms_comp = [rms_comp; sqrt(sum(rms_part(:,k))/size(rms_part,1))];
end

    % State-wise RMS
rms_part = [];
for k = 1:size(stateError,1)
    rms_part = [rms_part; norm(stateError(k,:))^2];
end
rms_full = sqrt(sum(rms_part)/length(rms_part));

end

