%% Prelab 2 Question 5
%   - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Constants
H = 8.773; % K/in
T0 = 278.083; % K
L = 5.875; % in
alpha = 4.82e-5 / (0.0254)^2; % m^2/s -> in^2/s

tMax = 1000; % sec
nMax = 10; % int

% xVec = 0:0.01:L;
x = 4.875; % in, location of Th8
tVec = 1:1:tMax;
nVec = 0:1:nMax;
nVec = reshape(nVec,1,1,[]);

%% u(x,t)
lambda = @(n) ((2*n - 1)*pi)/(2*L);
b = @(n) (-1).^n .* (8*H*L)./((2.*n - 1).^2*pi^2);
fourier = @(x,t,n) b(n).*sin(lambda(n).*x).*exp(-lambda(n).^2*alpha.*t);
u = @(x,t,n) T0 + H*x + sum(fourier(x,t,n),3);

%% Calculate solution for k n terms
for k = 1:length(nVec)
    if nVec(k) == 0
        uSol{k} = (T0 + H*x).*ones(1,length(tVec));
    else
        uSol{k} = u(x,tVec,nVec(2:k));
    end
end

%% Extract t = 1 sec and t = 1000 sec for k fourier terms
for k = 1:length(uSol)
    u1(k) = uSol{k}(1);
    u1000(k) = uSol{k}(1000);
end

%% Plot Convergence
figure
hold on
title("Heat Equation Convergence at Th_8 with n Fourier Terms, Discrete Times")
plot(reshape(nVec,1,[]), u1)
plot(reshape(nVec,1,[]), u1000)
xlabel("Number of Fourier Terms")
ylabel("Temperature at Th8 [K]")
legend("t = 1 sec", "t = 1000 sec", 'Location','best')

Fo = alpha*[tVec(1),tVec(end)]/L^2

figure
hold on
title("Heat Equation Convergence at Th_8 with n Fourier Terms, Continuous Times")
for k = 1:length(nVec)
    label(k) = sprintf("n = %.0f", k-1);
    plot(tVec, uSol{k}, 'LineWidth', 2 - 0.1*k)
end
xlabel("Time (sec)")
ylabel("Temperature at Th8 [K]")
legend(label,'Location','best')


