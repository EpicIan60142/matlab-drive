
%% 
% When uploading to MATLAB Grader do not include clear
% MATLAB Grader checks for variables. Variable names should match
% Since we are copy pasting the code to matlab grader, filename isnt
% important, as long as the variable name matches


angle = 0:0.1:2*pi;
sinAngle = sin(angle);
cosAngle = cos(angle);

sinCosAngle = sinAngle + cosAngle;

sinCosAngle = sinCosAngle';  % transpose, so that they are aware grader wont accept a row vector

% show plot - not required for grader


