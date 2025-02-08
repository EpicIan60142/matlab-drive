%% ASEN 2003 Lab 3
% Section 302 Group 37: Ben Chapel, Alyxis Ellington, Ian Faber

clc; clear; close all;

r = 7.5; % 7.5 cm
l = 26; % 26 cm
d = 15.5; % 15.5 cm
w = 2; % 1 rad/sec
delta_r = 0.05; % 0.5mm uncertainty
delta_d = 0.05; % 0.5mm uncertainty
delta_l = 0.05; % 0.5mm uncertainty

revs = 6;

theta = (linspace(0, revs*2*pi, 100))';

[v_mod, beta] = LCSMODEL(r, d, l, theta, w);

figure
sgtitle("Locomotive Crank Shaft Simulation Across 6 Revolutions")

subplot(1,2,1)
hold on;
title("Shaft Y Velocity vs. theta")
plot(theta, v_mod(:,2))
xlabel("Theta (rad)")
ylabel("Shaft Y Velocity (cm/s)")

subplot(1,2,2)
hold on;
title("Shaft X Velocity vs. theta")
plot(theta, v_mod(:,1))
xlabel("Theta (rad)")
ylabel("Shaft X Velocity (cm/s)")

theta_exp = [];
w_exp = [];
v_exp = [];
time = [];
files = {'Test1_5pt5V' 'Test1_6pt5V' 'Test1_7pt5V' 'Test1_8pt5V' 'Test1_9pt5V' 'Test1_10pt5V'};

res = cell(1,6);
f_dat = figure;
f_err = figure;
f_per_err = figure;
for i = 1:6
    [th, w, v, t] = LCSDATA(files{i});
    [v_mod, beta] = LCSMODEL(r, d, l, th, w);
    
    sigma = l^2 - (d-r*sin(th)).^2;
    part_r = -w.*(sin(th) + cos(th).*(d - r*sin(th))./sqrt(sigma) - (l^2*r*sin(th).*cos(th))./sigma.^(3/2));
    part_d = -r*w.*l^2.*cos(th)./(sigma).^(3/2);
    part_l = r*w.*l.*cos(th).*(d-r*sin(th))./sigma.^(3/2);
    delta_v = sqrt((delta_r*part_r).^2 + (delta_d*part_d).^2 + (delta_l*part_l).^2);
    
    figure(f_dat)
    subplot(2,3,i)
    hold on
    plot(th, v)
    errorbar(th, v_mod(:,2),delta_v)
    hold off
    grid on
    xlabel('Theta (rad)')
    ylabel('Slide Speed (cm/s)')
    title(['Slide Speed vs Crank Angle at ' num2str(i+4.5) 'V'])
    legend('Real', 'Model')
    ylim([-210 210])
    
    figure(f_err)
    res{i} = v - v_mod(:,2);
    subplot(2,3,i)
    hold on
    errorbar(th, res{i},delta_v,'r')
    plot(th, res{i},'b')
    hold off
    grid on
    xlabel('Theta (rad)')
    ylabel('Residuals (cm/s)')
    title(['Residuals vs Crank Angle at ' num2str(i+4.5) 'V'])
    ylim([-25 25])
    
    figure(f_per_err)
    subplot(2,3,i)
    plot(th, res{i}./v)
    grid on
    xlabel('Theta (rad)')
    ylabel('Residuals (% of measured value)')
    title(['Residuals vs Crank Angle at ' num2str(i+4.5) 'V'])
    ylim([-17 15])
end

res_avg = cellfun(@mean,res);
res_std = cellfun(@std,res);

