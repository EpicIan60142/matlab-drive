%% Project01 Problem 3.2.2

%% Setup
clc; clear;

h = 0.01;
simTime = 0:h:50; % 5001 elements

A_01 = zeros(1,length(simTime));
A_01(1) = 750000;
A_02 = zeros(1,length(simTime));
A_02(1) = 750000;

t_0 = 0;
t_1 = 0;
t_2 = 0;

p1 = 4000;
p2 = 4500;

interest1 = 0;
interest2 = 0;

paidOff1 = false;
paidOff2 = false;


%% Case 1: p = $4000/month
for k = 1:5001 %100*(simTime)+1, used to ensure positive indices for matrix assignment
    if(A_01(k) <= 0 && paidOff1 == false)
        A_01(1) = 750000;
        t_1 = h*(k-1); %Convert k back into years
        interest1 = 12*p1*t_1 - A_01(1);
        fprintf("Case 1: \nYears to repay: %f \nInterest paid: %f \n\n", t_1, interest1);
        paidOff1 = true;
    else
        A_01(k+1) = A_01(k) + eulerFunc(h*(k-1),A_01(k),p1)*h;
    end
end

%% Case 2: p = $4500/month
for k = 1:5001
    if(A_02(k) <= 0 && paidOff2 == false)
        A_02(1) = 750000;
        t_2 = h*(k-1);
        interest2 = 12*p2*t_2 - A_02(1);
        fprintf("Case 2: \nYears to repay: %f \nInterest paid: %f \n\n", t_2, interest2);
        paidOff2 = true;
    else
        A_02(k+1) = A_02(k) + eulerFunc(h*(k-1),A_02(k),p2)*h;
    end
end

%% Plotting
figure
hold on; grid on;
plot(0:h:t_1,A_01(1,1:100*t_1+1),'r-');
plot(0:h:t_2,A_02(1,1:100*t_2+1),'b-');
title("Time to pay off an adjustable rate mortgage based on monthly payments");
xlabel("Time to pay off (years)");
ylabel("Amount owed (USD)");
legend("$4000/month","$4500/month");

