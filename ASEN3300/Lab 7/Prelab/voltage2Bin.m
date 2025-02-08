function [binDec, binBin] = voltage2Bin(minV, maxV, bits, voltages)
% voltage2Bin: Determines the bin number, in both decimal and binary, 
%              a given voltage signal would be placed in by an A/D 
%              converter

range = maxV - minV; 
binSize = range/(2^bits);

binDec = floor((voltages-minV)/binSize);
binBin = dec2bin(binDec);

end