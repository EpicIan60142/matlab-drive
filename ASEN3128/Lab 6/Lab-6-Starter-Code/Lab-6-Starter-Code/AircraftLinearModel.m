function [Alon, Blon, Alat, Blat] = AircraftLinearModel(trim_definition, trim_variables, aircraft_parameters)
%
% STUDENT COMPLETE


u0 = trim_definition(1);
h0 = trim_definition(2);

alpha0 = trim_variables(1);
de0 = trim_variables(2);
dt0 = trim_variables(3);

theta0 = alpha0;

ap = aircraft_parameters;
rho = stdatmo(h0);



%%%%%%%%%%%%%%%%%%%%%%
%%% Longitudinal
%%%%%%%%%%%%%%%%%%%%%%

%%%% Trim values
CW0 = ap.W/((1/2)*rho*u0^2*ap.S);
CL0 = CW0*cos(theta0);
CD0 = ap.CDmin + ap.K*(CL0-ap.CLmin)^2;
CT0 = CD0 + CW0*sin(theta0);


%%%% Nondimensional stabiulity derivatives in body coordinates

%%%%% This is provided since we never discussed propulsion - Prof. Frew
dTdu = dt0*ap.Cprop*ap.Sprop*(ap.kmotor-2*u0+dt0*(-2*ap.kmotor+2*u0));
CXu = dTdu/(.5*rho*u0*ap.S)-2*CT0;

CDu = 0;%CDM*Ma;    % Compressibility only.  Ignore aeroelasticity and other effects (dynamic pressure and thrust).
CLu = 0;%CLM*Ma;
Cmu = 0;%CmM*Ma;

CZu = [];

CZalpha =[];
CXalpha = [];

CZalphadot = [];

CZq = [];


% Longitudinal dimensional stability derivatives (from Etkin and Reid)
Xu = [];
Zu = [];
Mu = [];


Xw = [];
Zw = [];
Mw = [];

Xq = 0;
Zq = [];
Mq = [];

Xwdot = 0;
Zwdot = [];
Mwdot = [];


% Matrices
Alon = [];

Blon = zeros(6,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lateral
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lateral-directional dimensional stability derivatives
Yv = [];
Yp = [];
Yr = [];

Lv = [];
Lp = [];
Lr = [];

Nv = [];
Np = [];
Nr = [];

G = ap.Ix*ap.Iz-ap.Ixz^2;

G3=ap.Iz/G;
G4=ap.Ixz/G;
G8=ap.Ix/G;


Alat = [];

Blat = zeros(6,2);


