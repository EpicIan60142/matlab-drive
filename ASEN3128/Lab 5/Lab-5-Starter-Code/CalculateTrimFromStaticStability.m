% Eric W. Frew and Zachary N. Sunberg
% ASEN 3128
% Created: 10/15/20
% STUDENTS COMPLETE THIS FUNCTION

function [alpha_trim, elevator_trim] = CalculateTrimFromStaticStability(trim_definition, aircraft_parameters)
%
% Inputs:	trim_definition         = [V0; h0]
%           aircraft_parameters     = structure with A/C parameters
%
%
% Outputs:	alpha_trim
%           elevator_trim
%
%
% Methodology: Uses linearized force and moment balance to estimate elevator and aoa for trim

ap = aircraft_parameters;
Va_trim = trim_definition(1);
h_trim= trim_definition(2);
rho_trim = stdatmo(h_trim);

% Students complete function below
% Determine lift coefficient needed for trim.
CL_trim = ap.W / (0.5*rho_trim*Va_trim^2*ap.S);

% Solve system of equations for angle of attack and elevator angle.
delta = ap.CLalpha*ap.Cmde - ap.CLde*ap.Cmalpha;

alpha_trim = (ap.Cm0*ap.CLde + ap.Cmde*CL_trim)/delta;
elevator_trim = -(ap.Cm0*ap.CLalpha + ap.Cmalpha*CL_trim)/delta;

