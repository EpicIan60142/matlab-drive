clear; clc; close all;

w0 = 2001;
L = 27.25;

discretized = discretize_load_Ian(4, L, w0);

forces = discretized(:,1)
distances = discretized(:,2)

[Ay, Ma] = wall_reactions(discretized)

errorSingle = moment_error(discretized, L, w0)

maxSections = 30;
for p = 1:maxSections
    matrix = discretize_load_Ian(p, L, w0); 
    error(p) = moment_error(matrix, L, w0);
end

figure
hold on;
title("Error between discretized load and analytical beam moments")
plot(1:p, zeros(length(1:p)), 'k--');
plot(error, 'b-');
xlabel("Number of sections")
ylabel("Error (lb-ft)")
