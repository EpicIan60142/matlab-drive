function phaseIdx = PRNPhase(PRN)
% Function that determines the correct phase bits to use depending on the
% provided PRN number
%   Inputs:
%       - PRN: PRN number to use
%   Outputs:
%       - phaseIdx: Array of bit indices used in the phase selector for CA
%                   generation
%
%   By: Ian Faber, 09/03/2025
%

switch PRN
    case 1
        idx = [2 6];
    case 2
        idx = [3 7];
    case 3
        idx = [4 8];
    case 4
        idx = [5 9];
    case 5
        idx = [1 9];
    case 6
        idx = [2 10];
    case 7
        idx = [1 8];
    case 8
        idx = [2 9];
    case 9
        idx = [3 10];
    case 10
        idx = [2 3];
    case 11
        idx = [3 4];
    case 12
        idx = [5 6];
    case 13
        idx = [6 7];
    case 14
        idx = [7 8];
    case 15
        idx = [8 9];
    case 16
        idx = [9 10];
    case 17
        idx = [1 4];
    case 18
        idx = [2 5];
    case 19
        idx = [3 6];
    case 20
        idx = [4 7];
    case 21
        idx = [5 8];
    case 22
        idx = [6 9];
    case 23
        idx = [1 3];
    case 24
        idx = [4 6];
    case 25
        idx = [5 7];
    case 26
        idx = [6 8];
    case 27
        idx = [7 9];
    case 28
        idx = [8 10];
    case 29
        idx = [1 6];
    case 30
        idx = [2 7];
    case 31
        idx = [3 8];
    case 32
        idx = [4 9];
    case 33
        idx = [5 10];
    case 34
        idx = [4 10];
    case 35
        idx = [1 7];
    case 36
        idx = [2 8];
    case 37
        idx = [4 10];
    otherwise
        fprintf("\nInvalid PRN requested, defaulting to PRN 1\n");
        idx = [2 6];
end

phaseIdx = idx;

end