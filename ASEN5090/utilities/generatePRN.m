function [CA, G1, G2] = generatePRN(G1, G2, PRN, chips)
% Calculates the PRN sequence for a given PRN number, number of chips, and 
% G1/G2 initialization
%   Inputs:
%       - G1: Initial G1 shift register contents, nominally ones(10,1)
%       - G2: Initial G2 shift register contents, nominally ones(10,1)
%       - PRN: PRN number of the code to generate, from 1 to 37
%       - chips: Number of chips to calculate, from 1 to 1023
%   Outputs:
%       - CA: specified-chip C/A code
%       - G1: Output of G1 register after the 1023 chip code
%       - G2: Output of G2 register after the 1023 chip code
%
%   By: Ian Faber, 09/03/2025
%

% Initialize CA
CA = zeros(1, chips);

% Determine phase
phaseIdx = PRNPhase(PRN);

% Generate 1023 chip CA code
for k = 1:length(CA)
    % G1 register
    out1 = G1(end);
    newBit1 = xor(G1(3), G1(10));
    G1 = [newBit1, G1(1:9)];

    % G2 register
    out2 = G2(end);
    newBit2 = xor(G2(2), xor(G2(3), xor(G2(6), xor(G2(8), xor(G2(9), G2(10))))));
    G2i = xor(G2(phaseIdx(1)), G2(phaseIdx(2)));
    G2 = [newBit2, G2(1:9)];

    % CA code
    CA(k) = xor(out1, G2i);

end


end