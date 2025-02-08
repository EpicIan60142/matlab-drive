function outputData = RWData(filename)
% Function that reads in data from various MEMS experiment result files
%   Inputs: Filename of result file
%   Outputs: outputData: Structure with fields of time, gyroOutput, and
%   inputRate

data = readmatrix(filename);
time = data(:,1)/1000;
cmdTorque = data(:,2)/1000; % Convert from mNm to Nm
angV = data(:,3);
actCurrent = data(:,4);

% Moving average to smooth noisy data
cmdTorque = (cmdTorque(1:end-4) + cmdTorque(2:end-3) + cmdTorque(3:end-2) + cmdTorque(4:end-1) + cmdTorque(5:end))/5;

outputData = struct('time', time, 'cmdTorque', cmdTorque, 'angV', angV, 'actCurrent', actCurrent);

end