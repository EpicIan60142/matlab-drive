%Inputs: Current time, current state
%Outputs: Rate of change of the state
%Ex. if y s velocity, output is acceleration

function out = odeFun(t,y)
    %make sure output is same size
    out = zeros(size(y));

    state1 = y(1);
    state2 = y(2);
    
    out1 = state1 + t; %velocity
    out2 = state2 - t; %acceleration
    
    out(1) = out1;
    out(2) = out2;
end
