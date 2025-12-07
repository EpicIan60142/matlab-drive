function [PRIF, iono] = ionocorr(rho1, f1, rho2, f2)
% Function that calculates the iono-free pseudorange combination and
% ionospheric correction based on pseudorange measurments on two different
% carrier frequencies
%   Inputs:
%       - rho1: Pseudorange data on carrier frequency 1 in meters
%       - f1: Carrier frequency 1 in Hz
%       - rho2: Pseudorange data on carrier frequency 2 in meters
%       - f2: Carrier frequency 2 in Hz
%   Outputs:
%       - PRIF: iono-free pseudorange combination in meters
%       - iono: Ionospheric correction for each carrier frequency, arranged
%               as follows: iono = [iono_1, iono_2]
%
%   By: Ian Faber, 10/03/2025
%

    % Calculate iono-free pseudorange combination
PRIF = ((f1^2)*rho1 - (f2^2)*rho2)/(f1^2 - f2^2);

    % Calculate total electron content
TEC = ((f1^2*f2^2)/(40.3*(f1^2 - f2^2)))*(rho2 - rho1);

    % Calculate ionosphere delay for each channel
iono = [40.3*TEC/(f1^2), 40.3*TEC/(f2^2)];


end

