function outputData = GyroData(filename)
% Function that reads in data from various MEMS experiment result files
%   Inputs: Filename of result file
%   Outputs: outputData: Structure with fields of time, gyroOutput, and
%   inputRate

data = readmatrix(filename);
time = data(2:end,1) - data(2,1);
gyroOutput = data(2:end, 2);
inputRate = data(2:end, 3);

outputData = struct('time', time, 'gyroOutput', gyroOutput, 'inputRate', inputRate);

end