% lab 5 nonlinear equations 


%% Constants

%syms zeta omegan Kp Kd

Kp = 1; % proprotial gain  3
Kd = 1; % dervitive gain -1.3
Kg = 33.3; % gear ratio
Km = 0.0401; % poprtional motor constant relating the speed to the motor voltage 
J  = 0.0005 + 0.2*(0.2794^2) + 0.0015; % moment of intertia 
Rm = 19.2; % output resitance of the motor 

%var = [Kp,Kd];

omegan = sqrt((Kp*Kg*Km)/(J*Rm));
ts   = 1;
zeta = (Kg^2 * Km^2  + Kd*Kg*Km)/(2*sqrt(Kp*Kg*Km*J*Rm));
%zeta = vpasolve(0.2  == exp(-((zeta*pi)/sqrt(1 - zeta^2))))
%omegan = vpasolve(0.05 == 1 -(1/sqrt(1- zeta^2))*exp(-zeta*omegan*ts)*sin(omegan*sqrt(1 - zeta^2)*ts + acos(zeta)))
%Kp = vpasolve(abs(omegan) == sqrt((Kp*Kg*Km)/(J*Rm)))
% Y = vpasolve(eqn,var)

% vpasolve([0.2  == exp(-((zeta*pi)/sqrt(1 - zeta^2))),...
  %   0.05 == 1 -(1/sqrt(1- zeta^2))*exp(-zeta*omegan*ts)*sin(omegan*sqrt(1 - zeta^2)*ts + acos(zeta))],[zeta,omegan]);
%0.05 == 1 -(1/sqrt(1- zeta^2))*exp(-zeta*omegan*ts)*sin(omegan*sqrt(1 - zeta^2)*ts + acos(zeta));
w=0;
for i = -10:0.5:10   % Kp values
    for h = -5:0.5:5 % Kd values
      if i == 0 
          break 
      else 
        w = w+1;
        Kp = i;
        Kd = h;
        omegan = sqrt((Kp*Kg*Km)/(J*Rm));
        ts   = 1;
        zeta = (Kg^2 * Km^2  + Kd*Kg*Km)/(2*sqrt(Kp*Kg*Km*J*Rm));
        xt  = 1 -(1/sqrt(1- zeta^2))*exp(-zeta*omegan*ts)*sin(omegan*sqrt(1 - zeta^2)*ts + acos(zeta));
        finalm(w,1) = i;
        finalm(w,2) = h;
        finalm(w,3) = xt;
      end
    end
end



