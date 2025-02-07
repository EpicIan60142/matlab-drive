%% Prelab 2 Question 6
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
nVec = 1;
alphaVec = (2e-5:1e-5:8e-5)/(0.0254)^2;
% alphaVec = reshape(alphaVec,1,1,[]);

%% u(x,t)
lambda = @(n) ((2*n - 1)*pi)/(2*L);
b = @(n) (-1).^n .* (8*H*L)./((2.*n - 1).^2*pi^2);
fourier = @(x,t,n,alpha) b(n).*sin(lambda(n).*x).*exp(-lambda(n).^2.*alpha.*t);
u = @(x,t,n,alpha) T0 + H*x + sum(fourier(x,t,n,alpha),3);

%% Calculate solution for k n terms
for k = 1:length(nVec)
    uSol{k} = u(x,tVec,nVec(1:k),alpha);
end

%% Extract t = 1 sec and t = 1000 sec for k fourier terms
for k = 1:length(uSol)
    u1(k) = uSol{k}(1);
    u1000(k) = uSol{k}(1000);
end

%% Plot Convergence
figure
hold on
title("Heat Equation Convergence at Th_8 with 1 Fourier Term")
for k = 1:length(nVec)
    plot(tVec, uSol{k}, 'LineWidth', 2 - 0.1*k)
end
xlabel("Time (sec)")
ylabel("Temperature at Th8 [K]")

%% Alpha trade study
for k = 1:length(alphaVec)
    uSol{k} = u(x,tVec,nVec,alphaVec(k));
end

figure
hold on
title("Heat Equation Convergence at Th_8 with 1 Fourier Term")
for k = 1:length(alphaVec)
    label(k) = sprintf("\\alpha = %.2d in^2/s", alphaVec(k));
    plot(tVec, uSol{k}, 'LineWidth', 2 - 0.1*k)
end
xlabel("Time (sec)")
ylabel("Temperature at Th8 [K]")
legend(label, 'Location','best')

