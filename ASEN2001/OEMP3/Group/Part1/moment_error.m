function error = moment_error(matrix, L, w0)
% Function that calculates the error between the discretized moment
% from discretized_load and the exact solution I found by hand
%   Inputs: Matrix of resultant forces and applied distances, length
%           of the beam L, and the force constant w0
%   Outputs: Error between the discretized moment and the analytical
%            moment

forces = matrix(:,1);
distances = matrix(:,2);

[~, Ma] = wall_reactions(matrix);

d = (3*L)/16;

M_dist = @(x)(w0/2)*(((-x^3)/(3*L))+(x^2)+(L*x)-((L^2)/3));
M_shear = @(x)(w0*x)*((L/2) + x - ((x^2)/(2*L)));

M_forces = 0;
for i = 1:length(distances)
    if(distances(i) < d)
        M_forces = M_forces + (forces(i)*distances(i));
    end
end

M_point = M_shear(d) - M_forces - Ma;

error = M_point - M_dist(d);

end

