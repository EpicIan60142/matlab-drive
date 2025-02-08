clc; clear; close all;

syms x
frac = 0.05;
tFinal = 8;
Kp = [7, 15, 16, 17, 64];
Kd = 0:2:10;
s = tf('s');

p = zeros(2, length(Kp));

for k = 1:length(Kp)
    % Part b
    G(k) = Kp(k) / (s^2 + 8*s + Kp(k));
    val = Kp(k);
    f = x^2 + 8*x + val;
    p(:,k) = solve(f, x);

    % Part c
    omega_n(k) = sqrt(Kp(k));
    zeta(k) = 4/omega_n(k);
    if zeta(k) < 1
        Mp(k) = exp(-(zeta(k)/sqrt(1 - zeta(k)^2)*pi));
        ts(k) = -log(frac*sqrt(1-zeta(k)^2))/(zeta(k)*omega_n(k));
    else
        Mp(k) = 0;
        ts(k) = 0;
    end

    figure
    hold on
    grid on
    titleText = sprintf("Unit Step Response of \\Theta(s) with K_p = %.0f", Kp(k));
    title(titleText);
    [y,t] = step(G(k));
    response = plot(t, y);
    xlabel("Time (sec)")
    ylabel("\theta (rad)")

    if zeta(k) < 1
        overshoot = yline(1+Mp(k), 'r');
        settleTime = xline(ts(k), 'b--');
        lowerP = yline(1-frac, 'k--');
        upperP = yline(1+frac, 'k--');
        
        lowerLabel = sprintf("%.1f%% Settling Time Bound", 100*frac);
        MpLabel = sprintf("M_p = %.3e rad", Mp(k));
        TsLabel = sprintf("t_s = %.3f sec", ts(k));

        legend("Response", MpLabel, TsLabel, lowerLabel, 'Location', 'best');
    end
    
end

Kp = [Kp, 1000];

figure
hold on
grid on
title("Unit Step Response of \Theta(s) with varying K_p, K_d = 0")
for k = 1:length(Kp)
    G(k) = Kp(k) / (s^2 + 8*s + Kp(k));
    [y, t] = step(G(k), tFinal);
    plot(t, y);
    label(k) = sprintf("K_p = %.0f", Kp(k));
end
xlabel("Time (sec)")
ylabel("\theta (rad)")

legend(label, 'Location', 'best')

% Part e/f
Kp0 = 64;
Kd = [Kd, 30];
Tfinal = 1.5;
figure
hold on
grid on
title("Unit Step Response of \Theta(s) with varying K_d, K_p = 64")
for k = 1:length(Kd)
    G2(k) = (Kd(k)*s + Kp0)/(s^2 + (Kd(k)+8)*s + Kp0);
    [y, t] = step(G2(k), Tfinal);
    plot(t,y)
    label(k) = sprintf("K_d = %.0f", Kd(k));
end
xlabel("Time (sec)")
ylabel("\theta (rad)")

legend(label, 'Location', 'best')
