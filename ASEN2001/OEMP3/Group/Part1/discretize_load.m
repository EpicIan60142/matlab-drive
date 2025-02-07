function [matrix] = discretize_load(p, L, w0)
% Function that splits up a distributed load into p sections,
% calculating the resultant force and location of each section
%   Inputs: Number of sections (p), length of the beam (L),
%           force constant w0
%   Outputs: p x 2 matrix of resultant forces and locations
%
% For this function, I am using the trapezoidal rule since it is
% the exact area for this situation and generally a good option
% for approximating. I am calculating xBar with the 1-D centroid 
% formula, i.e. integral(x*f(x)dx)/integral(f(x)dx), both integrals
% are calculated from x0 to xf, the start and stop of each section.
% I am assuming that the force distribution is acting according to
% a constant function, and does not change mathematically from section
% to section. If that assumption was not true, I would have another
% input that can accomodate different function handles.

sectionSize = L/p;

f = @(x) w0*(1-(x/L));

matrix = [];

for i = 1:p
    x0 = sectionSize*(i-1);
    xf = sectionSize*i;
    x = linspace(x0, xf, 1000);
    
    Fr(i) = trapz(x, f(x));
    xBar(i) = trapz(x, x.*f(x))/Fr(i);
    matrix = [matrix; [Fr(i), xBar(i)]];
end


end

