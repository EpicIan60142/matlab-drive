function [slope, intercept] = LeastSquares(xVector,yVector)
%Description: Calculates the slope and intercept of a best fit line for a
%specific number of x and y data points
%Inputs: x data point vector, y data point vector
%Outputs: Line slope, line intercept

A = sum(xVector); %Sum up the values of the x vector
B = sum(yVector); %Sum up the values of the y vector
C = sum(xVector.*yVector); %Sum up the product of x and y values from their specific vectors
D = sum(xVector.^2); %Sum up the squares of the values of the x vector
N = length(xVector); %Find the length of the x vector

slope = ((A*B)-(N*C))/((A^2)-(N*D)); %Calculate the slope of best fit
intercept = ((A*C)-(B*D))/((A^2)-(N*D)); %Calculate the intercept of best fit

end

