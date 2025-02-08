% Homework 9: Angle of Attack and Lift
clc %Clear the command window

%% Data Generation Section

angleVector = linspace(-5,25,100); %Generate a vector of 100 angles of attack between -5 and 25 degrees
yData = GenerateData(angleVector); %Call the GenerateData function to create synthetic experimental data 
figure %Create a new figure
plot(angleVector,yData,'b.'); %Plot the experimental data

%% Least Squares Section

index = 1; %Set the initial vector index to 1
while(index <= length(angleVector)) %Loop until all values of the angle vector are read
    if(angleVector(index) < 10) %If the value of the angle vector at an index is less than 10 degrees
        index1 = index; %Set the end of the first group of values to the current index
    elseif(angleVector(index) >= 10 && angleVector(index) < 15) %If the value of the angle vector at an index is between 10 and 15 degrees
        index2 = index; %Set the end of the second group of values to the current index
    elseif(angleVector(index) >= 15 && angleVector(index) < 20) %If the value of the angle vector at an index is between 15 and 20 degrees
        index3 = index; %Set the end of the third group of values to the current index
    else %If the value of the angle vector at an index is greater than 20 degrees
        index4 = index; %Set the end of the fourth and final group of values to the current index
    end
    
    index = index + 1; %Increment the index by 1
end

[M(1), B(1)] = LeastSquares(angleVector(1:index1),yData(1:index1)); %Call the LeastSquares function to calculate the first slope and intercept values
[M(2), B(2)] = LeastSquares(angleVector(index1+1:index2),yData(index1+1:index2)); %Call the LeastSquares function to calculate the second slope and intercept values
[M(3), B(3)] = LeastSquares(angleVector(index2+1:index3),yData(index2+1:index3)); %Call the LeastSquares function to calculate the third slope and intercept values
[M(4), B(4)] = LeastSquares(angleVector(index3+1:index4),yData(index3+1:index4)); %Call the LeastSquares function to calculate the fourth slope and intercept values

%% Piecewise Section

xVector = linspace(-5,25,1000); %Create a vector of 1000 angles of attack between -5 and 25 degrees

pindex = 1; %Set the initial plot index to 1
while(pindex <= length(xVector)) %Loop until all values of the x vector are read
    if(xVector(pindex) < 10) %If the angle of attack is less than 10 degrees
        pindex1 = pindex; %Set the first end index to the current index
    elseif(xVector(pindex) >= 10 && xVector(pindex) < 15) %If the angle of attack is between 10 and 15 degrees
        pindex2 = pindex; %Set the second end index to the current index
    elseif(xVector(pindex) >= 15 && xVector(pindex) < 20) %If the angle of attack is between 15 and 20 degrees
        pindex3 = pindex; %Set the third end index to the current index
    else %If the andle of attack is greater than 20 degrees
        pindex4 = pindex;%Set the final end index to the current index
    end
    
    pindex = pindex + 1; %Increment the plot index by 1
end

yFit = PiecewiseLeastFit(M,B,xVector); %Call the PiecewiseLeastFit function to generate y values for the line of bast fit
hold on %Save all plots
grid on %Turn on the coordinate grid
plot(xVector(1:pindex1),yFit(1:pindex1)); %Plot the first line of best fit
plot(xVector(pindex1+1:pindex2),yFit(pindex1+1:pindex2)); %Plot the second line of best fit
plot(xVector(pindex2+1:pindex3),yFit(pindex2+1:pindex3)); %Plot the third line of best fit
plot(xVector(pindex3+1:pindex4),yFit(pindex3+1:pindex4)); %Plot the fourth line of best fit
