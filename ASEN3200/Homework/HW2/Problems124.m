clc; clear; close all;

%% Problem 1
mass = [10, 10, 8, 8, 12, 12];
xM = [1, -1, 4, -2, 3, -3];
yM = [1, -1, -4, 2, -3, 3];
zM = [1, -1, 4, -2, -3, 3];

M = sum(mass);

xG = sum(mass.*xM)/M;
yG = sum(mass.*yM)/M;
zG = sum(mass.*zM)/M;

ITotal = zeros(3,3);

for k = 1:length(mass)
    x = xM(k) - xG;
    y = yM(k) - yG;
    z = zM(k) - zG;

    I{k} = [
                mass(k)*(y^2 + z^2), -mass(k)*x*y, -mass(k)*x*z;
                -mass(k)*x*y, mass(k)*(x^2 + z^2), -mass(k)*y*z;
                -mass(k)*x*z, -mass(k)*y*z, mass(k)*(y^2 + y^2)
           ];

%     disp(I{k});

    ITotal = ITotal + I{k};
end

disp(ITotal)

%% Problem 2
V = [1; 2; 2];
u = V/norm(V);

Iv = u'*ITotal*u;

%% Problem 4

[eigVectors, eigValues] = eig(ITotal);

FB = [
        eigVectors(:,1)';
        eigVectors(:,2)';
        eigVectors(:,3)'
     ];

IFrame = FB*ITotal*FB'












