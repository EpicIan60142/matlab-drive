%Josh Mellin
%This function checks if the sample rate is within acceptable bounds and
%returns the appropriate sample rate if it is not.  It takes the current
%value as well as the maximum allowable value.  Negative numbers aren't
%allowed, as well as decimal numbers.

function sr = check_sample_format(val, mval)
sr = val;
if(val > mval)
    sr = mval;
end
if(val < 1) %make sure it's not negative or less than 1
    sr = 1;
end
if(isnan(val)) %make sure they put a number in
    sr = 10;
end

sr = floor(sr); %take the floor in case someone put in a decimal

end


