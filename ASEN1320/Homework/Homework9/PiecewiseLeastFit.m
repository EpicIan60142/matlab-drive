function piecewiseVector = PiecewiseLeastFit(slopeVector, interceptVector, xVector)
%Description: Evaluates a linear piecewise function
%Inputs: Slopes vector, intercepts vector, x data points vector
%Outputs: Piecewise linear vector

piecewiseVector = zeros(1,length(xVector)); %Define a vector that is the same length as the input x vector

index = 1; %Initialize the index to 1
while(index <= length(xVector)) %Loop until all values of the x vector are read
    if(xVector(index) < 10) %If the xVector is less than 10
        piecewiseVector(index) = slopeVector(1)*xVector(index) + interceptVector(1); %Calculate y values based on the first line of best fit
    elseif(xVector(index) >= 10 && xVector(index) < 15) %If the xVector is between 10 and 15
        piecewiseVector(index) = slopeVector(2)*xVector(index) + interceptVector(2); %Calculate y values based on the second line of best fit
    elseif(xVector(index) >= 15 && xVector(index) < 20) %If the xVector is between 15 and 20
        piecewiseVector(index) = slopeVector(3)*xVector(index) + interceptVector(3); %Calculate y values based on the third line of best fit
    else %If the xVector is greater than 20
        piecewiseVector(index) = slopeVector(4)*xVector(index) + interceptVector(4); %Calculate y values based on the fourth line of best fit
    end
    
    index = index + 1; %Increment the index by 1
end

end

