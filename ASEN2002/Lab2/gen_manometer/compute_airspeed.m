%Josh Mellin

%This function computes the actual airspeed based on the measured
%pressures.  It's kind of the inverse of the compute_goal function.  It
%also returns the voltage based on the airspeed

function [airspeed, voltage] = compute_airspeed(pressure, ref1, ref2)
p = pressure(end);

%fill in the equations here
if(strcmp(ref1,'T1'))
    if(strcmp(ref2,'T2'))
        airspeed = p / 40;
    else if(strcmp(ref2,'T3'))
            airspeed = p/50;
        else if(strcmp(ref2,'T4'))
                airspeed = p/60;
            else
                airspeed = 0;
            end
        end
    end
end

if(strcmp(ref1,'T2'))
    if(strcmp(ref2,'T1'))
        airspeed = p / 40;
    else if(strcmp(ref2,'T3'))
            airspeed = p/50;
        else if(strmp(ref2,'T4'))
                airspeed = p/60;
            else
                airspeed = 0;
            end
        end
    end
end

if(strcmp(ref1,'T3'))
    if(strcmp(ref2,'T2'))
        airspeed = p / 40;
    else if(strcmp(ref2,'T1'))
            airspeed = p/50;
        else if(strmp(ref2,'T4'))
                airspeed = p/60;
            else
                airspeed = 0;
            end
        end
    end
end

if(strcmp(ref1,'T4'))
    if(strcmp(ref2,'T2'))
        airspeed = p / 40;
    else if(strcmp(ref2,'T3'))
            airspeed = p/50;
        else if(strmp(ref2,'T1'))
                airspeed = p/60;
            else
                airspeed = 0;
            end
        end
    end
end

voltage = airspeed / 100;





end







