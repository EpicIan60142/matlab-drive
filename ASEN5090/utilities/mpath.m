function [MP1, CMC1] = mpath(rho1, phi1, f1, phi2, f2)
% Function that calculates multi path observable and code minus carrier for
% a given signal
%   Inputs:
%       - rho1: Pseudorange for the signal of interest
%       - phi1: Carrier phase for the signal of interest
%       - f1: Frequency of the carrier signal of interest
%       - phi2: Carrier phase for the support signal
%       - f2: Frequency of the support carrier signal
%   Outputs:
%       - MP1: Multipath observable for the signal of interest
%       - CMC1: Code minus carrier for the signal of interest
%
%   By: Ian Faber, 10/03/2025
%

    % Speed of light
c = 299792458; % m/s

    % Calculate wavelength of carrier signals
lambda1 = c/f1;
lambda2 = c/f2;

    % Calculate multipath observable
MP1 = rho1 - ((f1^2 + f2^2)/(f1^2 - f2^2))*phi1*lambda1 + ((2*f2^2)/(f1^2 - f2^2))*phi2*lambda2;

    % Calculate code minus carrier
CMC1 = rho1 - phi1*lambda1;

end