function [Fc,Gc] = RotationDerivativeFeedback(initial_state_variables,m,g)
% Section 011 - Ian Faber, Felix Evrard, Blake Hogen, Gabriel Agostine
%         Function that calculates the control forces and moments needed
%         for derivative feedback control
% INPUTS: initial_state_variables is a vector of state variables at some
%                                 given time
%         m is the mass of the drone
%         g is acceleration due to gravity

k = 0.004;
Fc = [0;0;-m.*g];

L_c = -k.*initial_state_variables(10);
M_c = -k.*initial_state_variables(11);
N_c = -k.*initial_state_variables(12);
Gc = [L_c;M_c;N_c];

end