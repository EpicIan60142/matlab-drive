function [Ay, Ma] = wall_reactions(matrix)
% Function that calculates wall reactions based on a discretized
% distributed load
%   Inputs: Matrix of resultant forces and applied distances after
%           discretizing a distributed load
%   Outputs: Vertical wall reaction (Ay) and moment (Ma)

forces = matrix(:,1);
distances = matrix(:,2);

Ay = sum(forces);
Ma = sum(distances.*forces);

end

