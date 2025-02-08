clear; clc; close all;

RL = 10e6;
C = 100e-6;
I0 = 1;
tau = RL*C;

duty = 50;

current1 = @(t, width)I0*square((2*pi*t)/width, duty) > 0; % 8 tau
current2 = @(t, width)I0*square((2*pi*t)/width, duty) > 0; % tau
current3 = @(t, width)I0*square((2*pi*t)/width, duty) > 0; % tau / 2

const1.R = RL;
const1.C = C;
const1.width = 8*tau;
const1.current = current1;

const2.R = RL;
const2.C = C; 
const2.width = tau;
const2.current = current2;

const3.R = RL;
const3.C = C; 
const3.width = tau/2;
const3.current = current3;

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

X01 = 0;
X02 = 0;
X03 = 0;
tspan = 0:tau/1000:10*tau;

[t1, Vm1] = ode45(@(t,V)VoltageEOS(t,V,const1), tspan, X01, options);
[t2, Vm2] = ode45(@(t,V)VoltageEOS(t,V,const2), tspan, X02, options);
[t3, Vm3] = ode45(@(t,V)VoltageEOS(t,V,const3), tspan, X03, options);

time = 0:10*tau;
Vm = @(t)RL*I0 - RL*I0*exp(-t/tau);

figure
subplot(3,1,1)
hold on;
title("V_m(t) vs. time, Current pulse width of 8\tau")
yyaxis left
plot(t1, Vm1, 'r');
plot(time, Vm(time), 'k--');
ylabel("V_m(t) (V)");
yyaxis right
plot(time, current1(time, const1.width),'b-.');
ylim([0 10])
ylabel("Input current (A)");
xlabel("time (sec)");
legend("Actual membrane voltage", "Constant current membrane voltage", "Input current",'Location','best')


subplot(3,1,2)
hold on;
title("V_m(t) vs. Current, pulse width of \tau")
plot(t2, Vm2, 'm');
plot(time, Vm(time), 'k--');
ylabel("V_m(t) (V)");
yyaxis right
plot(time, current2(time, const2.width),'b-.');
ylim([0 10])
ylabel("Input current (A)");
xlabel("time (sec)");
legend("Actual membrane voltage", "Constant current membrane voltage", "Input current",'Location','best')


subplot(3,1,3)
hold on;
title("V_m(t) vs. Current, pulse width of \tau/2")
plot(t3, Vm3, 'g');
plot(time, Vm(time), 'k--');
ylabel("V_m(t) (V)");
yyaxis right
plot(time, current3(time, const3.width),'b-.');
ylim([0 10])
ylabel("Input current (A)");
xlabel("time (sec)");
legend("Actual membrane voltage", "Constant current membrane voltage", "Input current",'Location','best')

