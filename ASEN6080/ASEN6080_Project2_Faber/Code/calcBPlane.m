function [BdotR, BdotT, sig_R, sig_T, sig_RT, X_crossing, P_Bplane, STR2ECI, XPhi_BPlane, t_BPlane] = calcBPlane(XPhi_3SOI, t_3SOI, P_3SOI, pConst, DynFunc, opt)
% Function that calculates the expected Bplane crossing based on when a S/C
% passes the 3 SOI mark of a given planet.
%   Inputs:
%       - XPhi_3SOI: The state at the point where the spacecraft crosses 3
%                    SOI, organized as follows:
%                    [X; Y; Z; Xdot; Ydot; Zdot; ...; Phi_col], where
%                    Phi_col is a columnized STM
%       - t_3SOI: The time since the initial epoch at the point where the
%                 spacecraft crosses 3 SOI
%       - P_3SOI: State covariance at the point where the spacecraft
%                 crosses 3 SOI
%       - pConst: Planetary constants structure, as defined in
%                 getPlanetConst.m
%       - DynFunc: Function handle of the EOM to apply for this flyby as a
%                  function of t and X, where X includes an STM
%       - opt: ODE45 options to apply to integration
%   Outputs:
%       - BdotR: Component of the B plane crossing in the Rhat direction of
%                the Bplane frame
%       - BdotT: Component of the B plane crossing in the That direction of
%                the Bplane frame
%       - X_crossing: State at the B plane crossing in ECI coordinates
%       - P_Bplane: Propagated state covariance at the B plane in STR
%                   coordinates
%       - STR2ECI: DCM that converts STR coordinates to ECI coordinates
%       - XPhi_BPlane: Trajectory during the Bplane calculation
%       - t_BPlane: Time during the Bplane calculation
%
%   By: Ian Faber, 04/14/2025
%

    % Extract size of state
n = size(P_3SOI, 1);

    % Extract state
rVec = XPhi_3SOI(1:3);
vVec = XPhi_3SOI(4:6);

    % Calculate hyperbolic orbit parameters
r_inf = norm(rVec);
v_inf = norm(vVec);

        % Eccentricity
eVec = (1/pConst.mu_Earth)*((v_inf^2 - (pConst.mu_Earth/r_inf))*rVec - dot(rVec, vVec)*vVec);
e = norm(eVec);

        % Phat vector
Phat = eVec/e;

        % Specific angular momentum
hVec = cross(rVec, vVec);
h = norm(hVec);

        % What vector
What = hVec/h;

        % Semimajor axis
a = -pConst.mu_Earth/(v_inf^2 - (2*pConst.mu_Earth/r_inf));

        % Semilatus rectum
p = (h^2)/pConst.mu_Earth;

        % Semiminor axis
b = abs(a)*sqrt(e^2 - 1);

        % Shat vector
Shat = vVec/v_inf; % In ECI coords!

        % Nhat vector (aligned with Earth's north pole)
Nhat = [0; 0; 1];

        % That vector
That = cross(Shat, Nhat)/norm(cross(Shat, Nhat)); % In ECI coords!

        % Rhat vector
Rhat = cross(Shat, That); % In ECI coords!

        % B vector
% B = b*cross(Shat, What);
B = rVec - dot(rVec, vVec)*vVec;

        % DCM from STR to ECI frame
STR2ECI = [Shat, That, Rhat];

        % True anomaly
cNu = dot(rVec/r_inf, Phat); % cos(nu)

        % Hyperbolic anomaly
f = acosh(1 + (v_inf^2/pConst.mu_Earth)*((a*(1-e^2))/(1 + e*cNu)));

    % Calculate LTOF
LTOF = (pConst.mu_Earth/(v_inf^3))*(sinh(f) - f);

    % Integrate EOM to B-plane crossing
tspan = [t_3SOI, t_3SOI + LTOF];
XPhi_0 = [XPhi_3SOI(1:n); reshape(eye(n),n^2,1)];
[t_BPlane, XPhi_BPlane] = ode45(@(t,X)DynFunc(t,X), tspan, XPhi_0, opt);

    % Extract state and STM at Bplane crossing
X_crossing = XPhi_BPlane(end,1:n);
Phi_crossing = reshape(XPhi_BPlane(end,n+1:end), n, n);

    % Propagate covariance to B plane
P_Bplane = Phi_crossing*P_3SOI*Phi_crossing'; % In ECI coords!

    % Rotate covariance to STR frame
blkRot = blkdiag(STR2ECI', STR2ECI', 0);
P_Bplane = blkRot*P_Bplane*blkRot';

    % Assign outputs
BdotR = dot(B,Rhat);
BdotT = dot(B,That);

sig_R = sqrt(P_Bplane(3,3));
sig_T = sqrt(P_Bplane(2,2));
sig_RT = P_Bplane(2,3);

end