function dV = VoltageEOS1(t,V,const)

%Inputs: time vector t, state vector X [V; I]
%Outputs: derivative vector dV [dV/dt, dI/dt]

RL = const.R;
C = const.C;
width = const.width;
I = const.current;

dV = (I(t, width) - V/RL)/C;

end