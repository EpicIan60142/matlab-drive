function [score, cost, error, Mdist, Mpoint] = bestTree(p, matrix, w0, L)

M_dist = @(x)(w0/2)*(((-x^3)/(3*L))+(x^2)+(L*x)-((L^2)/3));
M_shear = @(x)(w0*x)*((L/2) + x - ((x^2)/(2*L)));

forces = matrix(:,1);
distances = matrix(:,2);
[~, Ma] = wall_reactions(matrix);

x = linspace(0,L,100);
cost = 500*(p^2);

error = 0;

Mdist = zeros(1,length(x));
Mpoint = zeros(1,length(x));

for i = 1:length(x)
    d = x(i);
    M_forces = 0;
    
    for j = 1:length(distances)
        if(distances(j) < d)
            M_forces = M_forces + (forces(j)*distances(j));
        end
    end
    
    M_point = M_shear(d) - M_forces - Ma;
    
    error = error + (M_dist(d) - M_point)^2;
    
    Mdist(i) = M_dist(d);
    Mpoint(i) = M_point;
end

error = (1/100)*error;

score = cost + error;

end

