function dChiVec = sigPointEOM(t, ChiVec, EOMFunc)
% Function that simultaneously propagates the EOM for all sigma points in
% a UKF for Stat OD problems
%   Inputs:
%       - t: Current integration time
%       - ChiVec: L*(2L + 1) sigma point vector from a UKF, organized as
%                 follows:
%                 [Chi_1; Chi_2; ...; Chi_2L+1], where Chi_i is an Lx1
%                 sigma point vector
%       - EOMFunc: Function handle for the equations of motion to be
%                  applied as a function of t and the state. Any other
%                  constants should be passed in when the handle is
%                  created, i.e. mu, J2, and J3.
%                  Ex. EOMFunc = @(t,X)orbitEOM_MuJ2(t,X,pConst.mu, ...
%                                                    pConst.J2, pConst.Ri)
%   Outputs:
%       - dChiVec: Rate of change vector for the provided sigma points
%                  [dChi_1; dChi_2; ...; dChi_2L+1]
%
%   By: Ian Faber, 03/15/2025
%

    % Find L, the size of the individual sigma point vectors
total = size(ChiVec,1);
L = round(max(roots([2, 1, -total]))); % total = 2*L^2 + L -> 2*L^2 + L - total = 0

    % Split up ChiVec into its original sigma point vectors
Chi = reshape(ChiVec, L, 2*L+1);

    % Loop through each vector and stack their respective rates of change
rates = [];
for k = 1:(2*L + 1)
    rates = [rates; EOMFunc(t, Chi(:,k))];
end

    % Assign output
dChiVec = rates;

end