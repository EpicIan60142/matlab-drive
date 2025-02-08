function [out] = bubbleSort(vector)
%bubbleSort - sort the vector in descending order!
%   input: unsorted vector
%   output: sorted vector

arrayLength = size(vector);

for i = 1:arrayLength
    for i = 1:(arrayLength-1)
        while vector(i) < vector(i+1)
            swap = vector(i);
            vector(i) = vector(i+1);
            vector(i+1) = swap;
        end        
    end
end

out = vector;

end
