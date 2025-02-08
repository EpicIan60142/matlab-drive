function [sigBR, omegBR] = calcError(sigBN, omegBN, RN, omegRN)
% Calculates the attitude and angular velocity tracking errors between two
% frames B and R.
%
%   Inputs:
%       - t: Time to compute error at
%       - sigBN: Current MRP set describing the nano-satellite's
%                orientation at time t
%       - omegBN: Current angular velocity vector of the nano-satellite at
%                 time t, in body coordinates
%       - RN: Current reference frame orientation DCM
%       - omegRN: Current angular velocity vector of the reference frame at
%                 time t, in inertial coordinates
%   Outputs:
%       - sigBR: MRP set describing the attitude tracking error of B
%                relative to R
%       - omegBR: angular velocity vector describing the angular velocity
%                 tracking error of B relative to R, in B coordinates

BN = MRP2DCM(sigBN);
sigBR = DCM2MRP(BN*RN', 1);

omegBR = BN*(BN'*omegBN - omegRN);


end