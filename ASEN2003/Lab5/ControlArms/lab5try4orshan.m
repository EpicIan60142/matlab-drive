clc; clear; close all;

A = readmatrix('April72022\April72022\rigid_5_kp_20_kd_0');
B = readmatrix('April72022\April72022\rigid_5_kp_15_kd_1pt5');
C = readmatrix('April72022\April72022\rigid_5_kp_10_kd_0');
D = readmatrix('April72022\April72022\rigid_5_kp_5_kd_0');

diffA = A(1,1);
diffA2 = A(1,2);
for  i = 1:(length(A))
    A(i,1) = (A(i,1) - diffA)/1000;
    A(i,2) = A(i,2) - diffA2;
end
startA = 3150;

diffB = B(1,1);
diffB2 = B(1,2);
for  i = 1:(length(B))
    B(i,1) = (B(i,1) - diffB)/1000;
    B(i,2) = B(i,2) - diffB2;
end 
startB = 2289;

diffC = C(1,1);
diffC2 = C(1,2);
for  i = 1:(length(C))
    C(i,1) = (C(i,1) - diffC)/1000;
    C(i,2) = C(i,2) - diffC2;
end 
startC = 2956;

diffD = D(1,1);
diffD2 = D(1,2);
for  i = 1:(length(D))
    D(i,1) = (D(i,1) - diffD)/1000;
    D(i,2) = D(i,2) - diffD2;
end
startD = 1408;


[E1,E2] = Aaron(20,0,  A(length(A)/2,2));
[F1,F2] = Aaron(15,1.5,B(length(B)/2,2));
[G1,G2] = Aaron(10,0,  C(length(C)/2,2));
[H1,H2] = Aaron( 5,0,  D(length(D)/2,2));

exp = [A(:,1),A(:,2),B(:,1),B(:,2),C(:,1),C(:,2),D(:,1),D(:,2)];
%[x2,t2] = impulse(sysTF);

%figure(1)
%plot(t1,x1)
%hold on 
%plot(t2,x2)
%title('Bode plot for time vs angular position')
%xlabel('Time (sec)')
%ylabel('Angular Position (rad)')
%legend('step response', 'Impulse resposne')


figure(1)
plot(exp(:,1)-exp(startA,1),exp(:,2))
hold on 
plot(E2,E1)
title('Bode plot for time vs angular position with Kp = 20 and Kd = 0')
xlabel('Time (sec)')
ylabel('Angular Position (rad)')
legend('Experimental results', 'Model results')
xlim([0,E2(end)])
hold off

figure(2)
plot(exp(:,3)-exp(startB,3),exp(:,4))
hold on 
plot(F2,F1)
title('Bode plot for time vs angular position with Kp = 15 and Kd = 1.5')
xlabel('Time (sec)')
ylabel('Angular Position (rad)')
legend('Experimental results', 'Model results')
xlim([0,F2(end)])
hold off

figure(3)
plot(exp(:,5)-exp(startC,5),exp(:,6))
hold on 
plot(G2,G1)
title('Bode plot for time vs angular position with Kp = 10 and Kd = 0')
xlabel('Time (sec)')
ylabel('Angular Position (rad)')
legend('Experimental results', 'Model results')
xlim([0,G2(end,:)])
hold off

figure(4)
plot(exp(:,7)-exp(startD,7),exp(:,8))
hold on 
plot(H2,H1)
title('Bode plot for time vs angular position with Kp = 5 and Kd = 0')
xlabel('Time (sec)')
ylabel('Angular Position (rad)')
legend('Experimental results', 'Model results')
xlim([0,H2(end)])
hold off 





%% lab 5 nonlinear equations 

% zeta = (Kp*Kg*Km)/(J*Rm)
% ts   = 1
% omegn = (Kg^2 * Km^2  + Kd*Kg*Km)/(2*sqrt(Kp*Kg*Km*J*Rm))
% 0.2  = exp ^-((3*pi)/sqrt(1 - zeta^2))
% 0.05 = 1 -(1/sqrt(1- zeta^2))*exp^(-zeta*omegn*ts) *sin(omegn*sqrt(1 - zeta^2)*ts + acos(zeta)) 
% Y = vpasolve(eqn,var)

fprintf('done')

function [x,t] = Aaron(Kp,Kd, step)

%% Constants
Kg = 33.3; % gear ratio
Km = 0.0401; % poprtional motor constant relating the speed to the motor voltage 
J  = 0.0005 + 0.2*(0.2794^2) + 0.0015; % moment of intertia 
Rm = 19.2; % output resistance of the motor 
 

%% Closed Loop System
d2 = 1;                                        % term connected to s^2
d1 = ((Kg^2 * Km^2)/(J*Rm)) + ((Kd*Kg*Km)/(J*Rm)); % term connected to s
d0 = (Kp*Kg*Km)/(J*Rm);                        % term connected to 1

num = (Kp*Kg*Km)/(J*Rm);
den = [d2 d1 d0];
sysTF = tf(num,den);

%% Step Response
tInput = 0:1/1000:5; 
stepInput = @(t,c) c*ones(length(t),1);
[x,t] = lsim(sysTF, stepInput(tInput, step),tInput);
end 

