function Rn = calcPRNAutoCorr(CA, nShift)
% Function that calculates the autocorrelation of a given PRN C/A code
% subject to a specified chip shift
%   Inputs:
%       - CA: C/A code to calculate autocorrelation for
%       - nShift: Number of chips to shift PRN C/A code by
%   Outputs:
%       - Rn: Autocorrelation of PRN C/A code subject to a shift of nShift
%
%   By: Ian Faber, 09/04/2025
%

% Shift PRN code
CA_unshift = CA;
if nShift ~= 0
    CA_shift = shiftCA(CA, nShift);
else
    CA_shift = CA_unshift;
end

% Convert from [0, 1] to [1, -1]
%CA_unshift = convPRNZeroOne2PosNeg(CA_unshift);
%CA_shift = convPRNZeroOne2PosNeg(CA_shift);

% Calculate Rn
Rn = (1/length(CA_shift))*sum(CA_unshift.*CA_shift);

end