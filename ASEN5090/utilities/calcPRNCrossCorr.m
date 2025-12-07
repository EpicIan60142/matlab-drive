function Rn = calcPRNCrossCorr(CA_k, CA_l, nShift)
% Function that calculates cyclic cross correlation between two C/A codes
% for a specified chip shift
%   Inputs:
%       - CA_k: First C/A code to compare
%       - CA_l: Second C/A code to compare
%       - nShift: Number of chips to shift CA_l by
%   Outputs:
%       - Rn: Cyclic cross correlation between CA_k and CA_l shifted by
%             nShift chips
%
%   By: Ian Faber, 09/04/2025
%

% Shift CA_l
if nShift ~= 0
    CA_l = shiftCA(CA_l, nShift);
end


% Convert from [0, 1] to [1, -1]
%CA_k = convPRNZeroOne2PosNeg(CA_k);
%CA_l = convPRNZeroOne2PosNeg(CA_l);

% Calculate Rn
Rn = (1/length(CA_l))*sum(CA_k.*CA_l);



end