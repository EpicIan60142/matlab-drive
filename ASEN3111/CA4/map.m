function y = map(x, x_min, x_max, out_min, out_max)
% map - Maps one set of numbers to another set of numbers. Can convert
%       scalars to vectors/matrices, but not vice versa. Instead, it will
%       map each element in the matrix/vector.
%
%   Inputs: 
%       x: Number to convert (scalar, vector, or matrix)
%       x_min: Minimum of the input
%       x_max: Maximum of the input
%       out_min: Desired minimum of the output
%       out_max: Desired maximum of the output
%
%   Outputs: 
%       y: Mapped output number
%
%   Example use:
%       To remap a number, for example 50, currently in the range of 
%       [0, 255] to a new range [0, 10], call the function like so:
%
%       y = map(50, 0, 255, 0, 10)
%
%       The output y will then be 1.9608.
%
%   Author: Ian Faber, 4/16/2023
%
%   Adapted from standard Arduino Library function "map()"
%

    y = ((x - x_min).*(out_max - out_min)) ./ (x_max - x_min) + out_min;

end