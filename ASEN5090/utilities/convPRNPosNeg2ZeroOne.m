function CA_conv = convPRNPosNeg2ZeroOne(CA)
% Function that converts a C/A code from [1, -1] representation to [0, 1]
% representation
%   Inputs:
%       - CA: CA code expressed in [1, -1] representation
%   Outputs:
%       - CA_conv: CA code expressed in [0, 1] representation
%
%   By: Ian Faber, 09/04/2025
%

CA_conv = (CA == 1)*0 + (CA == -1)*1;

end