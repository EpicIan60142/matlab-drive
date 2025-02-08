function experimentalData = GenerateData(inputData)
%Description: Generates synthetic wind tunnel experimental data based on
%given input angles of attack
%Inputs: Vector of input angles of attack
%Outputs: Row vector with noisy experimental data
rng(uint64(now*1000)); %Seed the random number generator using the current time

polyCoeff = load('Pcoef.dat', '-ascii'); %Read the "Pcoef.dat" file for polynomial coefficients

experimentalData = polyval(polyCoeff,inputData); %Generate experimental data values using the Pcoef coefficients and an input data vector

for index = 1:length(experimentalData) %Repeat for the length of the experimentalData vector
    experimentalData(index) = experimentalData(index) + (-0.5 + 1*rand(1,1)); %Add some random noise to the experimental data
end

end

