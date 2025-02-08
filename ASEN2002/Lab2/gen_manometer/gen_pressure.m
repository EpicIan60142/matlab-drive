%Josh Mellin

%This function creates the next iteration of the pressure vector.  It caps
%the ramp up at a certain rate

%Inputs: current pressure vector, desired airspeed, desired voltage, T1,
%T2, T3, T4, pause time, number of points

function new_pressure = gen_pressure(pressure, da, ref1, ref2, pt, np, off, atm_press, rho)

last_point = pressure(end);
goal = compute_goal(ref1, ref2, da, atm_press, rho)+off;
thresh = 500; %maximum pressure differential to be achieved in 1 second
difference = last_point-goal;

if(abs(difference) < thresh*pt)
    new_pressure = linspace(pressure(end),goal,np) + randn(np,1)';
else
    if(difference > 0)
        new_pressure = linspace(pressure(end),pressure(end)-thresh*pt,np) + randn(np,1)';
    else
        new_pressure = linspace(pressure(end),pressure(end)+thresh*pt,np) + randn(np,1)';
    end
end







end









