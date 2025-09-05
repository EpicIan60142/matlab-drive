function CA_conv = convPRNZeroOne2PosNeg(CA)
% Function that converts a C/A code from [0, 1] representation to [1, -1]
% representation
%   Inputs:
%       - CA: CA code expressed in [0, 1] representation
%   Outputs:
%       - CA_conv: CA code expressed in [1, -1] representation
%
%   By: Ian Faber, 09/04/2025
%

CA_conv = (CA == 0)*1 + (CA == 1)*-1;

end