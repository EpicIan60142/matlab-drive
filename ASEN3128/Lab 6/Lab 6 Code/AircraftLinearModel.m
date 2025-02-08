function [Alon, Blon, Alat, Blat] = AircraftLinearModel(trim_definition, trim_variables, aircraft_parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:	trim_definition         = [V0; h0]
%           trim_variables          = [alpha0; de0; dt0]
%           aircraft_parameters     = structure with A/C parameters
%
%
% Outputs:	Alon [6x6]
%           Blon [6x2]
%           Alat [6x6]
%           Blat [6x2]
%
%
% Develop Longitudinal and Lateral A and B matrices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

CZu = 0;

CZalpha = -CD0 - ap.CLalpha;
CXalpha = CL0*(1-2*ap.K*ap.CLalpha);

CZalphadot = -ap.CLalphadot;

CZq = -ap.CLq;


% Longitudinal dimensional stability derivatives (from Etkin and Reid)
Xu = rho*u0*ap.S*CW0*sin(theta0) + 0.5*rho*u0*ap.S*CXu;
Zu = -rho*u0*ap.S*CW0*cos(theta0) + 0.5*rho*u0*ap.S*CZu;
Mu = 0.5*rho*u0*ap.S*ap.c*Cmu;


Xw = 0.5*rho*u0*ap.S*CXalpha;
Zw = 0.5*rho*u0*ap.S*CZalpha;
Mw = 0.5*rho*u0*ap.S*ap.c*ap.Cmalpha;

Xq = 0;
Zq = 0.25*rho*u0*ap.c*ap.S*CZq;
Mq = 0.25*rho*u0*ap.c^2*ap.S*ap.Cmq;

Xwdot = 0;
Zwdot = 0.25*rho*ap.c*ap.S*CZalphadot;
Mwdot = 0.25*rho*ap.c^2*ap.S*ap.Cmalphadot;


% Matrices
m = ap.m;
g = ap.g;
Iy = ap.Iy;

Alon = [
            Xu/m,                              Xw/m,                               0,                                          -g*cos(theta0),                           0,  0;
            Zu/(m-Zwdot),                      Zw/(m-Zwdot),                      (Zq + m*u0)/(m-Zwdot),                       (-m*g*sin(theta0))/(m-Zwdot),             0,  0;
            (1/Iy)*(Mu+(Mwdot*Zu)/(m-Zwdot)), (1/Iy)*(Mw + (Mwdot*Zw)/(m-Zwdot)), (1/Iy)*(Mq + (Mwdot*(Zq + m*u0)/(m-Zwdot))), (-Mwdot*m*g*sin(theta0))/(Iy*(m-Zwdot)),  0,  0;
            0,                                 0,                                  1,                                          0,                                        0,  0;
            1,                                 0,                                  0,                                          0,                                        0,  0;
            0,                                 1,                                  0,                                          -u0,                                      0,  0
       ];

Blon = zeros(6,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lateral
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lateral-directional dimensional stability derivatives
Yv = 0.5*rho*u0*ap.S*ap.CYbeta;
Yp = 0.25*rho*u0*ap.b*ap.S*ap.CYp;
Yr = 0.25*rho*u0*ap.b*ap.S*ap.CYr;

Lv = 0.5*rho*u0*ap.S*ap.b*ap.Clbeta;
Lp = 0.25*rho*u0*ap.S*ap.b^2*ap.Clp;
Lr = 0.25*rho*u0*ap.S*ap.b^2*ap.Clr;

Nv = 0.5*rho*u0*ap.S*ap.b*ap.Cnbeta;
Np = 0.25*rho*u0*ap.S*ap.b^2*ap.Cnp;
Nr = 0.25*rho*u0*ap.S*ap.b^2*ap.Cnr;

G = ap.Ix*ap.Iz-ap.Ixz^2;

G3=ap.Iz/G;
G4=ap.Ixz/G;
G8=ap.Ix/G;




Alat = [Yv/m                    Yp/m                  (Yr/m)-u0          ap.g*cos(theta0)       0                 0;
        G3*Lv + G4*Nv           G3*Lp + G4*Np         G3*Lr + G4*Nr      0                      0                 0;
        G4*Lv + G8*Nv           G4*Lp + G8*Np         G4*Lr + G8*Nr      0                      0                 0;
        0                       1                     tan(theta0)        0                      0                 0;
        0                       0                     sec(theta0)        0                      0                 0;
        1                       0                     0                  0                      u0*cos(theta0)    0];

Blat = zeros(6,2);


