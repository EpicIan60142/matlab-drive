function [dX, fGrav, fDrag, Fthrust, Pair] = rocketEOM(t, IC, const)
x=IC(1);
z=IC(2);
vx=IC(3);
vz=IC(4);
m=IC(5);
mair=IC(6);
vair=IC(7);

if (m < const.m_bottle)
    m=const.m_bottle + mair;
end

Abottle = pi*(const.d_bottle/200)^2; %convert diameter from cm to m
Athroat = pi*(const.d_throat/200)^2; %convert diameter from cm to m

v=[vx; vz];

P_amb = const.P_amb*6894.76;
Vi_air = const.V_bottle-const.Vi_water;
Pi_air = (const.P_gage+const.P_amb)*6894.76; 
rhoi_air = Pi_air/(const.Ti_air*const.R); 
mi_air = rhoi_air*Vi_air;

rho_air = mair/vair;

if (norm([x (z-const.z0)]) <  const.l_stand)
    headingstate="ONSTAND";
else
    headingstate="INFLIGHT";
end

if vair < const.V_bottle && rho_air > const.rho_air && z > 0
    flightstate="WATERTHRUST";
elseif vair >= const.V_bottle && rho_air > const.rho_air && z > 0
    flightstate="AIRTHRUST";
    if vair > const.V_bottle
        vair=const.V_bottle;
    end
elseif vair >= const.V_bottle && rho_air <= const.rho_air && z > 0
    flightstate="BALLISTIC";
    if rho_air < const.rho_air
        rho_air=const.rho_air;
    end
else
    flightstate="GROUND";
end

switch headingstate
    case "ONSTAND"
        heading=[cosd(const.theta); sind(const.theta)];
    case "INFLIGHT"
        heading=v/norm(v);
    otherwise
        heading=[0; 0];
end

switch flightstate
    case "WATERTHRUST"
        fGrav=[0; -m*const.g];
        fDrag= -heading*(.5*const.C_drag*Abottle*rho_air*norm(v)^2);
        Pair=Pi_air*(Vi_air/vair)^const.gamma;
        dmdt=-const.Cd*Athroat*sqrt(2*const.rho_water*(Pair-P_amb));
        dmAdt=0;
        dVdt=const.Cd*Athroat*sqrt((2/const.rho_water)*(Pi_air*((Vi_air/vair)^const.gamma)-P_amb));
        Fthrust=heading*(2*const.Cd*Athroat*(Pair-P_amb));
    case "AIRTHRUST"
        fGrav=[0; -m*const.g]; %possibly wrong because not right weight
        fDrag=-heading*(.5*const.rho_air*const.C_drag*Abottle*norm(v)^2);
        Pair=Pi_air*((mair*Vi_air)/(mi_air*vair))^const.gamma;
        
        if (Pair < P_amb)
            Pair = P_amb;
        end
        
        Tair=Pair/(rho_air*const.R);
        Pcrit=Pair*((2/(const.gamma+1))^(const.gamma/(const.gamma-1)));
              
        if (Pcrit > P_amb)
            Mach=1;
            Te=(2/(const.gamma+1))*Tair;
            Pe=Pcrit;
            rhoE=Pe/(const.R*Te);
        else 
            Mach=sqrt((2/(const.gamma-1))*(((Pair/P_amb)^((const.gamma-1)/const.gamma))-1));
            Te=Tair/(1+((const.gamma-1)/2)*Mach^2);
            Pe=P_amb;
            rhoE=P_amb/(const.R*Te);
        end
        
        Ve=Mach*sqrt(const.gamma*const.R*Te);
        
        dmAdt=-const.Cd*rhoE*Athroat*Ve;
        dmdt=-const.Cd*rhoE*Athroat*Ve;
        dVdt=0;
        Fthrust=heading*(-dmAdt*Ve+Athroat*(P_amb-Pe));
        
    case "BALLISTIC"
        Pair=P_amb;
        fGrav=[0; -m*const.g];
        fDrag=-heading*(.5*const.C_drag*Abottle*const.rho_air*(norm(v)^2));
        Fthrust=[0;0];
        dmdt=0;
        dmAdt=0;
        dVdt=0;

    otherwise
        Pair=P_amb;
        vx=0;
        vz=0;
        fGrav=[0;0];
        fDrag=[0;0];
        Fthrust=[0;0];
        dmdt=0;
        dmAdt=0;
        dVdt=0;
        
end
fNet=Fthrust+fDrag+fGrav;
ax=fNet(1)/m;
az=fNet(2)/m;
dX=[vx; vz; ax; az; dmdt; dmAdt; dVdt];

fprintf("Time: %f Heading State: %s , Flight State: %s , Thrust: [%f, %f], V_Air: %f , Rho_Air: %f, Z: %f , Mass_Air: %f , Air Pressure: %f \n", t, headingstate, flightstate, Fthrust(1), Fthrust(2), vair, rho_air, z, mair, Pair)
end
