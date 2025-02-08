clc; clear; close all;
% [w_mod, mom_F] = UWMODEL(theta,w)

files = {'Data\balanced_1','Data\balanced_2','Data\unbalanced_1','Data\unbalanced_2'};

mom_F = 0;
mom_F2 = 3;

[time1, theta1, w1] = UWData(files{1});
[time2, theta2, w2] = UWData(files{2});
[time3, theta3, w3] = UWData(files{3});
[time4, theta4, w4] = UWData(files{4});

theta = (0:15/35:15)';
[w1_test, ~, ~] = UWMODEL(1,theta,w1,mom_F);
[~, w2_check, ~] = UWMODEL(2,theta,w1_test,mom_F2);
[~, w3_test, ~] = UWMODEL(3,theta,w1_test,mom_F2);
[~, w4_test, ~] = UWMODEL(4,theta,w1_test,mom_F2);

figure()
hold on
title("Expected \omega vs. \theta")
plot(theta, w1_test)
plot(theta, w2_check)
plot(theta, w3_test)
plot(theta, w4_test)
hold off
xlabel('\theta (rad)')
ylabel('\omega (rad/s)')
legend('model 1','model 2','model 3','model 4')

[w1_mod, ~, ~] = UWMODEL(1,theta1,w1,mom_F);
figure()
hold on 
plot(theta1,w1)
plot(theta1,w1_mod)
hold off
xlabel('theta (rad)')
ylabel('Omega (rad/s)')
legend('Real','Model')
grid on

[~, w2_mod, mom_F] = UWMODEL(2,theta2,w2,mom_F2);
figure()
hold on 
plot(theta2,w2)
plot(theta2,w2_mod)
hold off
xlabel('theta (rad)')
ylabel('Omega (rad/s)')
legend('Real','Model')
grid on

[~, w2_check1, mom_Ft1] = UWMODEL(2,theta2,w2,0);
[~, w2_check2, mom_Ft2] = UWMODEL(2,theta2,w2,-.5);
[~, w2_check3, mom_Ft3] = UWMODEL(2,theta2,w2,-1);
[~, w2_check4, mom_Ft4] = UWMODEL(2,theta2,w2,-1.5);
[~, w2_check5, mom_Ft5] = UWMODEL(2,theta2,w2,-2);
[~, w2_check6, mom_Ft6] = UWMODEL(2,theta2,w2,-2.5);
[~, w2_check7, mom_Ft7] = UWMODEL(2,theta2,w2,mom_Ft1);
figure()
hold on 
plot(theta2,w2)
plot(theta2,w2_check1)
plot(theta2,w2_check2)
plot(theta2,w2_check3)
plot(theta2,w2_check4)
plot(theta2,w2_check5)
plot(theta2,w2_check6)
plot(theta2+.15,w2_check7)
hold off
xlabel('theta (rad)')
ylabel('Omega (rad/s)')
legend('Real','0','-.5','-1','-1.5','-2','-2.5','calc')
grid on

mom_F = -1.5;
[~, w1_mod1, ~] = UWMODEL(1,theta1,w1,mom_F);
[~, w2_mod1, ~] = UWMODEL(1,theta2,w2,mom_F);
[~, w1_mod2, ~] = UWMODEL(2,theta1,w1,mom_F);
[~, w2_mod2, ~] = UWMODEL(2,theta2,w2,mom_F);
figure()
subplot(2,1,1)
hold on
plot(theta1,w1,'o')
plot(theta1,w1_mod1,'--')
plot(theta1,w1_mod2)
hold off
grid on
xlabel('theta (rad)')
ylabel('Omega (rad/s)')
legend('Real','Model 1','Model 2')
title('Balanced Test 1')

subplot(2,1,2)
hold on
plot(theta2,w2,'o')
plot(theta2,w2_mod1,'--')
plot(theta2,w2_mod2)
hold off
grid on
xlabel('theta (rad)')
ylabel('Omega (rad/s)')
legend('Real','Model 1','Model 2')
title('Balanced Test 2')

res1_mod1 = w1 - w1_mod1;
res1_mod2 = w1 - w1_mod2;
res2_mod1 = w2 - w2_mod1;
res2_mod2 = w2 - w2_mod2;
figure()
hold on
plot(theta1,res1_mod1)
plot(theta1,res1_mod2)
plot(theta2,res2_mod1)
plot(theta2,res2_mod2)
hold off
grid on
xlabel('theta (rad)')
ylabel('Omega Residuals (rad/s)')
legend('Test 1 Model 1','Test 1 Model 2','Test 2 Model 1','Test 2 Model 2')
title('Balanced Test Residuals')

res1_mod1_data = rescalc(res1_mod1)
res1_mod2_data = rescalc(res1_mod2)
res2_mod1_data = rescalc(res2_mod1)
res2_mod2_data = rescalc(res2_mod2)

[~, w3_mod3, ~] = UWMODEL(3,theta3,w3,mom_F);
[~, w4_mod3, ~] = UWMODEL(3,theta4,w4,mom_F);
[~, w3_mod4, ~] = UWMODEL(4,theta3,w3,mom_F);
[~, w4_mod4, ~] = UWMODEL(4,theta4,w4,mom_F);
figure()
subplot(2,1,1)
hold on
plot(theta3,w3,'o')
plot(theta3,w3_mod3)
plot(theta3,w3_mod4)
hold off
grid on
xlabel('theta (rad)')
ylabel('Omega (rad/s)')
legend('Real','Model 3','Model 4')
title('Unbalanced Test 1')

subplot(2,1,2)
hold on
plot(theta4,w4,'o')
plot(theta4,w4_mod3)
plot(theta4,w4_mod4)
hold off
grid on
xlabel('theta (rad)')
ylabel('Omega (rad/s)')
legend('Real','Model 3','Model 4')
title('Unbalanced Test 2')

res3_mod3 = w3 - w3_mod3;
res3_mod4 = w3 - w3_mod4;
res4_mod3 = w4 - w4_mod3;
res4_mod4 = w4 - w4_mod4;
figure()
hold on
plot(theta3,res3_mod3)
plot(theta3,res3_mod4)
plot(theta4,res4_mod3)
plot(theta4,res4_mod4)
hold off
grid on
xlabel('theta (rad)')
ylabel('Omega Residuals (rad/s)')
legend('Test 3 Model 3','Test 3 Model 4','Test 4 Model 3','Test 4 Model 4')
title('Unbalanced Test Residuals')

res3_mod3_data = rescalc(res3_mod3)
res3_mod4_data = rescalc(res3_mod4)
res4_mod3_data = rescalc(res4_mod3)
res4_mod4_data = rescalc(res4_mod4)

function calcs = rescalc(res)
    stdn = std(res);
    avg = mean(res);
    N = length(res);
    sdm = stdn/sqrt(N);
    Ng3s = sum(res > 3*stdn + avg | res < avg - 3*stdn);
    calcs = [stdn avg sdm N Ng3s];
end