clc; clear; close all;

A = readmatrix('Kp6_99KdN_121'); % 6.99,-0.121
B = readmatrix('Kp7KdN_5');   % 7, -0.5
C = readmatrix('Kp10KdN_5');   % 10, -0.5
D = readmatrix('Kp8Kd_01');     % 8,0.01

% DIFFS CORRESPOND TO
diffA = A(5480,1);
diffA2 = A(5480,2);
for  i = 1:(length(A))
    A(i,1) = (A(i,1) - diffA)/1000;
    A(i,2) = A(i,2) - diffA2;
end
startA = 5480;% hardcoded where data went to zero 

diffB = B(6566,1);
diffB2 = B(6566,2);
for  i = 1:(length(B))
    B(i,1) = (B(i,1) - diffB)/1000;
    B(i,2) = B(i,2) - diffB2;
end 
startB = 6566; % hardcoded where data went to zero 

diffC = C(5753,1);
diffC2 = C(5753,2);
for  i = 1:(length(C))
    C(i,1) = (C(i,1) - diffC)/1000;
    C(i,2) = C(i,2) - diffC2;
end 
startC = 5753; % hardcoded where data went to zero 

diffD = D(5198,1);
diffD2 = D(5198,2);
for  i = 1:(length(D))
    D(i,1) = (D(i,1) - diffD)/1000;
    D(i,2) = D(i,2) - diffD2;
end
startD = 5198; % hardcoded where data went to zero 


[E1,E2] = Aaron(6.99,-0.121, A(length(A),2));
[F1,F2] = Aaron(7, -0.5, B(length(B),2));
[G1,G2] = Aaron(10, -0.5, C(length(C),2));
[H1,H2] = Aaron(8,0.01, D(length(D),2));

%% section with for loops for Bs
for i = 0:5
h = 2*i+1;
g = 2*i+2;
[P1,P2] = Chris(6.99,-0.121,A(length(A),2),0.01*i);
Chrismat1(:,h)= P1;
Chrismat1(:,g)= P2;
end 

for i = 0:5
 h = 2*i+1;
 g = 2*i+2;
[P1,P2] = Chris(7,-0.5,B(length(B),2),0.01*i);
Chrismat2(:,h)= P1;
Chrismat2(:,g)= P2;
end 

for i = 0:5
h = 2*i+1;
    g = 2*i+2;
[P1,P2] = Chris(10,-0.5,C(length(C),2),0.01*i);
Chrismat3(:,h)= P1;
Chrismat3(:,g)= P2;
end 

for i = 0:5
    h = 2*i+1;
    g = 2*i+2;
[P1,P2] = Chris(8,0.01,D(length(D),2),0.01*i);
Chrismat4(:,h)= P1;
Chrismat4(:,g)= P2;
end 


%% Making larger data matrix
exp = [A(:,1),A(:,2),B(:,1),B(:,2),C(:,1),C(:,2),D(:,1),D(:,2)];

%% Plotting
figure(1)
plot(exp(:,1)-exp(startA,1),exp(:,2))
hold on 
plot(E2,E1)
title('Bode plot for time vs angular position with Kp = 6.99 and Kd = -0.121')
xlabel('Time (sec)')
ylabel('Angular Position (rad)')
legend('Experimental results', 'Model results','Location','best')
xlim([0,E2(end)])
hold off

figure(2)
plot(exp(:,3)-exp(startB,3),exp(:,4))
hold on 
plot(F2,F1)
title('Bode plot for time vs angular position with Kp = 7 and Kd = -0.5')
xlabel('Time (sec)')
ylabel('Angular Position (rad)')
legend('Experimental results', 'Model results','Location','best')
xlim([0,F2(end)])
hold off

figure(3)
plot(exp(:,5)-exp(startC,5),exp(:,6))
hold on 
plot(G2,G1)
title('Bode plot for time vs angular position with Kp = 10 and Kd = -0.5')
xlabel('Time (sec)')
ylabel('Angular Position (rad)')
legend('Experimental results', 'Model results','Location','best')
xlim([0,G2(end,:)])
hold off

figure(4)
plot(exp(:,7)-exp(startD,7),exp(:,8))
hold on 
plot(H2,H1)
title('Bode plot for time vs angular position with Kp = 8 and Kd = 1')
xlabel('Time (sec)')
ylabel('Angular Position (rad)')
legend('Experimental results', 'Model results','Location','best')
xlim([0,H2(end)])
hold off 

figure(5)%based on kp = 8 kd = 0.01
plot(exp(:,7)-exp(startD,7),exp(:,8))
hold on 
plot(Chrismat4(:,2),Chrismat4(:,1))
plot(Chrismat4(:,4),Chrismat4(:,3))
plot(Chrismat4(:,6),Chrismat4(:,5))
plot(Chrismat4(:,8),Chrismat4(:,7))
plot(Chrismat4(:,10),Chrismat4(:,9))
plot(Chrismat4(:,12),Chrismat4(:,11))
xlabel('Time (sec)')
ylabel('Angular Position (rad)')
legend('Exprimetal results','B=0','B=0.01','B=0.02','B=0.03','B=0.04','B=0.05','Location','best')
title('System Response for time vs angular position with Kp = 8 and Kd = 0.01')
xlim([0 4])
hold off

figure(6)%based on kp =6.99 kd = -0.121
plot(exp(:,1)-exp(startA,1),exp(:,2))
hold on 
plot(Chrismat1(:,2),Chrismat1(:,1))
plot(Chrismat1(:,4),Chrismat1(:,3))
plot(Chrismat1(:,6),Chrismat1(:,5))
plot(Chrismat1(:,8),Chrismat1(:,7))
plot(Chrismat1(:,10),Chrismat1(:,9))
plot(Chrismat1(:,12),Chrismat1(:,11))
xline(1,'k--');
yline(0.95*A(length(A),2),'r--')
yline(1.2*A(length(A),2),'b--')
yline(1.05*A(length(A),2),'r--')
xlabel('Time (sec)')
ylabel('Angular Position (rad)')
legend('Experimental results','B=0','B=0.01','B=0.02','B=0.03','B=0.04','B=0.05','Settling time of 1 sec','5% boundary of steady state value','20% overshoot boundary','Location','best')
title('System response for time vs angular position with Kp = 6.99 and Kd = -0.121')
xlim([0 4])
hold off

figure(7)%based on kp = 7 kd = -0.5
plot(exp(:,3)-exp(startB,3),exp(:,4))
hold on 
plot(Chrismat2(:,2),Chrismat2(:,1))
plot(Chrismat2(:,4),Chrismat2(:,3))
plot(Chrismat2(:,6),Chrismat2(:,5))
plot(Chrismat2(:,8),Chrismat2(:,7))
plot(Chrismat2(:,10),Chrismat2(:,9))
plot(Chrismat2(:,12),Chrismat2(:,11))
xline(1,'k--');
yline(0.95*B(length(B),2),'r--')
yline(1.2*B(length(B),2),'b--')
yline(1.05*B(length(B),2),'r--')
xlabel('Time (sec)')
ylabel('Angular Position (rad)')
legend('Experimental results','B=0','B=0.01','B=0.02','B=0.03','B=0.04','B=0.05','Settling time of 1 sec','5% boundary of steady state value','20% overshoot boundary','Location','best')
title('System response for time vs angular position with Kp = 7 and Kd = -0.5')
xlim([0 4])
hold off





%% Function creations

fprintf('done')

function [x,t] = Aaron(Kp,Kd, step)

%% Constants
Kg = 33.3; % gear ratio
Km = 0.0401; % poprtional motor constant relating the speed to the motor voltage 
J  = 0.0005 + 0.2*(0.2794^2) + 0.0015; % moment of intertia 
Rm = 19.2; % output resitance of the motor 
 

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

function [x,t] = Chris(Kp,Kd,step,B)

%% Constants
Kg = 33.3; % gear ratio
Km = 0.0401; % poprtional motor constant relating the speed to the motor voltage 
J  = 0.0005 + 0.2*(0.2794^2) + 0.0015; % moment of intertia 
Rm = 19.2; % output resitance of the motor 
 

%% Closed Loop System
d2 = 1;                                        % term connected to s^2
d1 = ((B/J)+(Kg^2 * Km^2)/(J*Rm)) + ((Kd*Kg*Km)/(J*Rm)); % term connected to s
d0 = (Kp*Kg*Km)/(J*Rm);                        % term connected to 1

num = (Kp*Kg*Km)/(J*Rm);
den = [d2 d1 d0];
sysTF = tf(num,den);

%% Step Response
tInput = 0:1/1000:5; 
stepInput = @(t,c) c*ones(length(t),1);
[x,t] = lsim(sysTF, stepInput(tInput, step),tInput);
end 
