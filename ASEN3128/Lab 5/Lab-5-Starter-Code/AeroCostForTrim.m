% Eric W. Frew
% ASEN 3128
% AeroCostForTrim.m
% Created: 10/15/20
% STUDENTS COMPLETE THIS FUNCTION

function cost = AeroCostForTrim(trim_variables, trim_definition, aircraft_parameters)
%
%
% INPUT:    trim_definition = [V0; h0]
%
%           trim_variables = [alpha0; de0; dt0];
%
% OUTPUT:   cost = norm(total_force) + norm(total_moment)
%
% 
% METHOD:   Determines the total force acting on the aircraft from the
%           aerodynamics and weight. Then takes the norm of both to create 
%           a single cost function that can be minimized.


rho0=stdatmo(trim_definition(2));

% Determine the TOTAL force `forces` and TOTAL moment `moments`
% acting on the aircraft based on the `trim_variables` and
% `trim definition` arguments. You should use the
% `AeroForcesAndMoments_BodyState_WindCoeffs` function

wind_inertial=[0,0,0]';
aircraft_state=[0,0,-trim_definition(2),0,trim_variables(1),0,trim_definition(1)*cos(trim_variables(1)),0,trim_definition(1)*sin(trim_variables(1)),0,0,0]';
aircraft_surfaces=[trim_variables(2),0,0,trim_variables(3)]';
[forces, moments, wind_angles] = AeroForcesAndMoments_BodyState_WindCoeffs(aircraft_state, aircraft_surfaces, wind_inertial, rho0, aircraft_parameters);

forces = forces + TransformFromInertialToBody([0; 0; aircraft_parameters.W], aircraft_state(4:6));

% Final cost is calculated from total force and moment vectors
cost = norm(forces) + norm(moments); 
end
