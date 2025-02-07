%Josh Mellin

%This function computes the goal pressure based on the desired airspeed and
%chosen pressure sensors.  It also takes the atmospheric pressure [pa] and
%air density [kg/m^3] as inputs

%P1 = Settling chamber static pressure
%P2 = Inlet static ring pressure
%P3 = Pitot-static probe total pressure port - the one at the front of the
%probe
%P4 = Pitot-static probe static pressure port - the ring around the probe

function goal = compute_goal(ref1, ref2, da, atm_press, rho)

%Equations for the pressures
if(strcmp(ref1, 'T1'))
    p1 = atm_press;
end

if(strcmp(ref1, 'T2'))
    p1 = atm_press - (1/2)*rho*da^2;
end

if(strcmp(ref1, 'T3'))
    p1 = atm_press;%*0.99;
end

if(strcmp(ref1, 'T4'))
    p1 = atm_press - (1/2)*rho*da^2;
end

if(strcmp(ref1, 'Vac'))
    p1 = 0;
end

if(strcmp(ref1, 'Amb'))
    p1 = atm_press;
end

%and set the second ref signal
if(strcmp(ref2, 'T1'))
    p2 = atm_press;
end

if(strcmp(ref2, 'T2'))
    p2 = atm_press - (1/2)*rho*da^2;
end

if(strcmp(ref2, 'T3'))
    p2 = atm_press;%*0.99;
end

if(strcmp(ref2, 'T4'))
    p2 = atm_press - (1/2)*rho*da^2;
end

if(strcmp(ref2, 'Vac'))
    p2 = 0;
end

if(strcmp(ref2, 'Amb'))
    p2 = atm_press;
end

%and then find the difference as the goal
goal = p1 - p2;

end







