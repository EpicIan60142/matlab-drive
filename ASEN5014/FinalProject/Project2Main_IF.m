%% ASEN 5014 Project 2 Main Script
% By: Ian Faber
% Partner: Gabriel Agostine

%% Housekeeping
clc; clear; close all;

%% Setup
format shortE

fprintf("--- Setup ---\n")

digits = 16; % Number of digits to keep for numerical precision

muEarth = 398600; % km^3/s^2
rOrbit = 6778; % km - 400 km LEO orbit

n = sqrt(muEarth/(rOrbit^3));
T = (2*pi)/n;

A = [
        0       0       0       1       0   0;
        0       0       0       0       1   0;
        0       0       0       0       0   1;
        3*n^2   0       0       0       2*n 0;
        0       0       0       -2*n    0   0;
        0       0       -n^2    0       0   0
    ];

B = [
        0 0 0;
        0 0 0;
        0 0 0;
        1 0 0;
        0 1 0;
        0 0 1
    ];

C = [
        1 0 0 0 0 0;
        0 1 0 0 0 0;
        0 0 1 0 0 0
    ];

D = [
        0 0 0;
        0 0 0;
        0 0 0
    ];

%% Modes
[G, Lambda] = eig(A); % Eigenvectors and eigenvalues

    % Fix numerical precision
G = round(G,digits);
Lambda = round(Lambda, digits);

if det(G) == 0 % G is singular, need more eigenvectors
    [V, J] = jordan(A); % Extended eigenvectors and jordan form
    X1 = V(:,1);
    X2 = V(:,2);
    U3 = real(V(:,3));
    V3 = imag(V(:,3));
    U4 = real(V(:,5));
    V4 = imag(V(:,5));
else
    X1 = G(:,1);
    X2 = G(:,2);
    U3 = real(G(:,3));
    V3 = imag(G(:,3));
    U4 = real(G(:,5));
    V4 = imag(G(:,5));
end

    % Eigenspaces and real modal spaces
eigenSpaces = V;
modeVecs = [X1, X2, U3, V3, U4, V4]

realModalForm = (modeVecs^-1)*A*modeVecs;

realModalForm = round(realModalForm, digits);

%% Closed loop design
    % Choose time constant
steadytime = linspace(600, 1500, 6);
tau = steadytime/5;

    % Create system poles
eigVals_cl = -1./tau

    % Design feedback controller
K_cl = place(A, B, eigVals_cl);

Acl = A-(B*K_cl);

F_cl = (C*((-A+(B*K_cl))^-1)*B)^-1;

    % xDot = (A-BK)x + BFr
    % y = (C-DK)x + DFr, D = zeros so y = Cx as before
B_F = B*F_cl;

%% Part 1
fprintf("--- Part 1 ---\n")
    % Observability Matrix
O = [C; C*A; C*A^2; C*A^3; C*A^4; C*A^5]

obsRank = rank(O)

    % Simulate unit perturbations in each modal space
t = 0:1:T;

phi = @(A, t) expm(A*t);
x = zeros(6,length(t),4);
y = zeros(3,length(t),4);
energy = zeros(4,1);

for k = 1:4
    switch k
        case 1
            x0 = X1/norm(X1);
        case 2
            x0 = X2/norm(X2);
        case 3
            x0 = U3/norm(U3);
        case 4
            x0 = U4/norm(U4);
        otherwise
            x0 = zeros(6,1); % Should never reach this
    end

        % Calculate energy of output in this modal space
    eng = [];
    for kk = 1:length(t)
        x(:,kk,k) = phi(A, t(kk))*x0;
        y(:,kk,k) = C*x(:,kk,k);
        eng = [eng; y(:,kk,k)'*y(:,kk,k)];
    end

    energy(k) = sum(eng);

end

    % Higher energy -> more observable
energy

%     % Find observability Grammian
% t01 = T;
% 
%     % Perturb A by a tiny amount for lyap
% Aperturb = A - 1e-8*eye(size(A)); %1e-10*diag([1, 1, 1, 0, 0, 0]);
% 
%     % Build Q
% Q = expm(-Aperturb*t01)*(C'*C)*expm(-Aperturb'*t01) - C'*C;
% G = lyap(Aperturb, Q)
% 
% obsEigVals = eig(G)

%% Part 2
fprintf("--- Part 2 ---\n")

    % a. Slow observer: observer eigenvalues = 0.2*closed loop eigenvalues
eigVals_Lslow = 0.2*eigVals_cl;
L_slow = place(A',C',eigVals_Lslow);
L_slow = L_slow'
eig_slow = eig(A-L_slow*C)

    % b. Equal observer: observer eigenvalues = closed loop eigenvalues
eigVals_Leq = eigVals_cl;
L_reg = place(A',C',eigVals_Leq);
L_reg = L_reg'
eig_reg = eig(A - L_reg*C)

    % c. Fast observer: observer eigenvalues = 5*closed loop eigenvalues
eigVals_Lfast = 5*eigVals_cl;
L_fast = place(A',C',eigVals_Lfast);
L_fast = L_fast'
eig_fast = eig(A - L_fast*C)

%% Part 3
fprintf("--- Part 3 ---\n")

    % Make observer cell array and set timescale
L = {L_slow, L_reg, L_fast};
t01 = T;

    % Step response options
stepOpt = timeoptions("cstprefs");
stepOpt.YLabel.String = "State Response";
stepOpt.Grid = 'on';
stepOpt.YLim = [0 1];
% stepOpt.OutputLabels.FontSize = 10;

    % Simulate unit step references
for k = 1:length(L)
    Aprime = [
                A-B*K_cl                                      B*K_cl
                zeros(size(A-L{k}*C,1), size(A-B*K_cl,2))     A-L{k}*C
             ];
    Bprime = [B_F; zeros(size(A-L{k}*C,1),size(B_F,2))];
    Cprime = [C, zeros(size(C,1), size(A-L{k}*C,1))];
    Dprime = zeros(size(C,1), size(B_F,2));
    sys = ss(Aprime, Bprime, Cprime, Dprime, 'StateName', {'x','y','z','vx','vy','vz', 'ex', 'ey', 'ez', 'evx', 'evy', 'evz'}, 'InputName', {'u1','u2','u3'}, 'OutputName', {'x [km]','y [km]','z [km]'});

    figure
    step(sys, t01, stepOpt)
    switch k
        case 1
            title("Closed loop slow observer/controller step response to unit reference inputs")
        case 2
            title("Closed loop regular observer/controller step response to unit reference inputs")
        case 3
            title("Closed loop fast observer/controller step response to unit reference inputs")
    end


end

%% Part 4
fprintf("--- Part 4 ---\n")

orthBasis = eye(12); %[IC_plant, 0; 0, IC_observer]

response = [];

for k = 1:length(L)
    Aprime = [
                A-B*K_cl                                      B*K_cl
                zeros(size(A-L{k}*C,1), size(A-B*K_cl,2))     A-L{k}*C
             ];

    responseComp = struct('L', L{k}, 'run', []);

    for kk = 1:size(orthBasis,1)
        x0 = orthBasis(:,kk); 

        x = zeros(size(Aprime,1), length(t));
        y = zeros(size(Cprime,1), length(t));
        for ii = 1:length(t)
            x(:,ii) = phi(Aprime, t(ii))*x0;
            y(:,ii) = Cprime*x(:,ii);
        end

        run = struct('x0', x0, 'x', x, 'y', y);

        responseComp.run = [responseComp.run; run];
    end

    response = [response; responseComp];

end

for k = 1:size(orthBasis,2)
    figure; tiles = tiledlayout(4,3);
    titleText = sprintf("System unit initial condition response,\nx_0 = [%.0f, %.0f, %.0f, %.0f, %.0f, %.0f, %.0f, %.0f, %.0f, %.0f, %.0f, %.0f]^T", orthBasis(:,k));
    title(tiles, titleText);
    for ii = 1:size(Aprime,1)
        nexttile
        hold on; grid on;  
        for kk = 1:length(L)
            switch kk
                case 1 % slow
                    slow = plot(t, response(kk).run(k).x(ii,:), 'b-');
                case 2 % regular
                    reg = plot(t, response(kk).run(k).x(ii,:), 'r-');
                case 3 % fast
                    fast = plot(t, response(kk).run(k).x(ii,:), 'm-');
            end
        end
        xlabel("Time [sec]")

        switch ii
            case 1
                ylabel("X [km]")
            case 2
                ylabel("Y [km]")
            case 3
                ylabel("Z [km]")
            case 4
                ylabel("X_{dot} [km/s]")
            case 5
                ylabel("Y_{dot} [km/s]")
            case 6
                ylabel("Z_{dot} [km/s]")
            case 7
                ylabel("e_X [km]")
            case 8
                ylabel("e_Y [km]")
            case 9
                ylabel("e_Z [km]")
            case 10
                ylabel("e_{X_{dot}} [km/s]")
            case 11
                ylabel("e_{Y_{dot}} [km/s]")
            case 12
                ylabel("e_{Z_{dot}} [km/s]")
        end

        if ii == 6
            legend([slow, reg, fast], ["Slow observer", "Regular observer", "Fast observer"], 'location', 'bestoutside')
        end

        
    end
    
end


%% Part 5
fprintf("--- Part 5 ---\n")

% For LQR:
%   Need (A,B) to be stabilizable -> all unstable modes reachable, (A,B) is
%   completely reachable -> (A,B) stabilizable

%   Need (A,V) to be detectable -> all unstable modes observable, (A,C) is
%   completely observable -> (A,C) detectable -> V = C, Q = V'*V = C'*C

Q = C'*C;
r = 5e8;
R = r*eye(3);

sys = ss(A,B,C,D);

K_LQR = lqr(sys, Q, R);


eigVals_LQR = eig(A-B*K_LQR)
tau_LQR = -1./real(eigVals_LQR')
tau

F_LQR = (C*((-A+(B*K_LQR))^-1)*B)^-1

B_F_LQR = B*F_LQR;

Aprime = [
            A-B*K_LQR                                       B*K_LQR
            zeros(size(A-L_fast*C,1), size(A-B*K_LQR,2))    A-L_fast*C
         ];
Bprime = [B_F_LQR; zeros(size(A-L_fast*C,1),size(B_F_LQR,2))];
Cprime = [C, zeros(size(C,1), size(A-L_fast*C,1))];
Dprime = zeros(size(C,1), size(B_F_LQR,2));

sys_LQR = ss(Aprime, Bprime, Cprime, Dprime, 'StateName', {'x','y','z','vx','vy','vz', 'ex', 'ey', 'ez', 'evx', 'evy', 'evz'}, 'InputName', {'u1','u2','u3'}, 'OutputName', {'x [km]','y [km]','z [km]'});

sys_cl = ss(A-B*K_cl, B_F, C, D, 'StateName', {'x','y','z','vx','vy','vz'}, 'InputName', {'u1','u2','u3'}, 'OutputName', {'x [km]','y [km]','z [km]'});

figure
stepplot(sys_cl, 5*max(tau), stepOpt)
title("System closed loop step response to unit inputs - pole placement controller")

[~,t_cl,x_cl] = step(sys_cl, 3*max(tau), stepOpt);

figure
stepOpt.YLim = [0 1.1];
stepplot(sys_LQR, 5*max(tau_LQR), stepOpt)
title("System closed loop step response to unit inputs - LQR")

[~,t_LQR,x_LQR] = step(sys_LQR, 3*max(tau_LQR), stepOpt);

ref = eye(3);

    % Pole placement energy
energy_pp = [];
u_cl = zeros(3,length(t_cl),3);
for kk = 1:3
    eng = [];
    for k = 1:length(t_cl)
        u = -K_cl*reshape(x_cl(k,:,kk)',6,1) + F_cl*ref(:,kk);
        eng = [eng; u'*u];
        u_cl(:,k,kk) = u;
    end
    eng = sum(eng);
    energy_pp = [energy_pp; eng];
end
energy_pp

    % LQR energy
u_LQR = zeros(3,length(t_LQR),3);
energy_LQR = [];
for kk = 1:3
    eng = [];
    for k = 1:length(t_LQR)
        u = -K_LQR*reshape(x_LQR(k,1:6,kk),6,1) + F_LQR*ref(:,kk);
        eng = [eng; u'*u];
        u_LQR(:,k,kk) = u;
    end
    eng = sum(eng);
    energy_LQR = [energy_LQR; eng];
end
energy_LQR

    % Plot control signals
for k = 1:3
    titleText = sprintf("Control signal comparison with r = [%.0f, %.0f, %.0f]^T \nover 3 time constants", ref(:,k));
    figure; tiles = tiledlayout(3,1);
    title(tiles, titleText)
    for kk = 1:3
        nexttile;
            hold on; grid on;
            titleText = sprintf("u_%.0f vs. time", kk);
            labelText = sprintf("u_%.0f [km/s^2]", kk);
            title(titleText)
            plot(t_cl, u_cl(kk,:,k));
            plot(t_LQR, u_LQR(kk,:,k));
            xlabel("Time [sec]"); ylabel(labelText)
    end
    legend("Pole Placement Control", "LQR Control", 'Location', 'bestoutside');
end



%% Reset format to default
format default;
