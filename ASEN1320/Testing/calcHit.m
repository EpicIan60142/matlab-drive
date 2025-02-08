function groundTime = calcHit(distanceVector,timeVector)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    i = 1;
    zeroFound = 0;
    while(i < length(timeVector) && zeroFound == 0)
        if(distanceVector(i) <= 0)
            groundTime = timeVector(i-1);
            zeroFound = 1;
        end
        i = i + 1;
    end
end

