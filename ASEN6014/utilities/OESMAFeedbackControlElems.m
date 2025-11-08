function dX = OESMAFeedbackControlElems(t, X, doe_r, kConst, pConst)
% Function that implements a semi-major axis only feedback control law for
% relative motion between a chief and deputy spacecraft, with desired 
% motion specified with orbit element differences.
%   Inputs:
%       - t: Current integration time in sec
%       - X: Current combined state vector as follows:
%            X = [X_c; X_d]
%       - doe_r: Structure of desired orbit element differences with any of 
%                the following fields, not necessarily all: 
%                   da, de, di, dRAAN, dargPeri, dM0
%                If a field is missing, it is ignored and assumed not
%                desirable to control.
%       - kConst: Structure of control gains with the following fields:
%                 - K: semi major axis gain
%       - pConst: Planetary constants vector containing mu for the
%                 celestial body of interest
%   Outputs:
%       - dX: Rate of change vector as follows:
%             [dX_c; dX_d]
%
%   By: Ian Faber, 11/06/2025

% Extract states
X_c = X(1:6);
X_d = X(7:12);

% Calculate state elements
cart.rVec = X_c(1:3);
cart.vVec = X_c(4:6);
cart.mu = pConst.mu;
oe_c = convCart2ClassicOE(cart);

cart.rVec = X_d(1:3);
cart.vVec = X_d(4:6);
cart.mu = pConst.mu;
oe_d = convCart2ClassicOE(cart);

rDeputy = cart.rVec;
vDeputy = cart.vVec;

% Calculate DCM
rHat = rDeputy/norm(rDeputy);
hHat = cross(rDeputy, vDeputy)/norm(cross(rDeputy, vDeputy));
thetaHat = cross(hHat, rHat);

NHd = [rHat, thetaHat, hHat];

% Calculate reference elements
oe_r.mu = pConst.mu;
oe_r.a = oe_c.a + doe_r.da;
oe_r.e = oe_d.e;
oe_r.i = oe_d.i;
oe_r.RAAN = oe_d.RAAN;
oe_r.argPeri = oe_d.argPeri;
oe_r.f = convE2f(convM2E(convf2E(oe_d.f, oe_d.e), oe_d.e, false), oe_r.e);

% Assign B matrix
cart = convClassicOE2Cart(oe_r);
r = norm(cart.rVec);
h = norm(cross(cart.rVec, cart.vVec));
p = oe_r.a*(1-oe_r.e^2);
eta = sqrt(1-oe_r.e^2);

B = [
        (2*oe_r.a*oe_r.e*sin(oe_r.f))/(h),              (2*oe_r.a^2*p)/(h*r),                0;
        (p*sin(oe_r.f))/h,                              ((p+r)*cos(oe_r.f) + r*oe_r.e)/h,    0;
        0,                                              0,                                   (r*cos(oe_r.argPeri + oe_r.f))/h;
        0                                               0,                                   (r*sin(oe_r.argPeri + oe_r.f))/(h*sin(oe_r.i));
        -(p*cos(oe_r.f))/h,                             ((p+r)*sin(oe_r.f))/(h*oe_r.e),      -(r*sin(oe_r.argPeri + oe_r.f)*cos(oe_r.i))/(h*sin(oe_r.i));
        (eta*(p*cos(oe_r.f)-2*r*oe_r.e))/(h*oe_r.e),    -(eta*(p+r)*sin(oe_r.f))/(h*oe_r.e), 0
    ];

% Assign gains
K = kConst.K;

% Calculate control effort and convert to inertial frame
da = oe_d.a - oe_r.a;
u = -K*B'*[da; zeros(5,1)];

u = NHd*u;

% Assign output
fChief = orbitEOM(t, X_c, pConst, struct(), false(1,2));
fChief = fChief(4:6);

fDeputy = orbitEOM(t, X_d, pConst, struct(), false(1,2));
fDeputy = fDeputy(4:6);

dX_c = [X_c(4:6); fChief];
dX_d = [X_d(4:6); fDeputy + u];

dX = [dX_c; dX_d];

end