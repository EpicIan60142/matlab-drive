function error = LSError(F, b, t)
% Function that calculates the least squares error of a regression line
% over some interval of the line
%   Inputs: Anonymous function handle F, data vector b, time interval t
%   Outputs: Least squares error

error = sqrt((sum((b-F(t)).^2))/(length(t)-2));

end
