function CA_shifted = shiftCA(CA, nShift)
% Function that shifts a C/A code by a specified number of chips
%   Inputs:
%       - CA: Unshifted C/A code
%       - nShift: number of chips to shift the C/A code by
%   Outputs:
%       - CA_shifted: C/A code shifted by nShift
%
%   By: Ian Faber, 09/04/2025
%

CA_unshift = CA;
CA_shifted = CA_unshift([((end-nShift+1):end) (1:end-nShift)]);

end