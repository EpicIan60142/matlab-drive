%% ASEN 2004 Lab 1 Milestone 2 Design Analyzer


%% Housekeeping
clc; clear; close all;

%% Data extraction

dataSet1 = 'Tempest UAS & B747 Airfoil and CFD Data for ASEN 2004 Aero Lab (Spr22).xlsx';
dataSet2 = 'Milestone2_2D_FlatPlate_Data.xlsx';

dataSet = dataSet2;

% Extract names and number of Excel sheets
[sheetStatus, sheetNames] = xlsfinfo(dataSet);
numSheets = length(sheetNames);

% Extract data
for k = 1:numSheets
    sheetData{k} = xlsread(dataSet, sheetNames{k});
end

if strcmp(dataSet, dataSet1)
    % Tempest UAS
    Data = cell2mat(sheetData(1,2));
    
    ClAlphas = Data(:,1);
    Cl = Data(:,2);
    Cd = Data(:,3);
    Re = Data(1,5);
else
    % Indexing variables (to match up Cl and Cd vectors, synced between -13 and 13 degrees AoA)
    startCl = 1;
    excludeCl = [5, 8, 10, 11, 13, 15, 16];
    stopCl = 24;
    startCd = 1;
    excludeCd = [1, 11, 20, 21, 22, 23];
    stopCd = 23;
    
    % Extract 2D coefficients
    ClData = cell2mat(sheetData(1));
    CdData = cell2mat(sheetData(2));
    
    ClAlphasRaw = ClData(startCl:stopCl,1);
    ClAlphas = ClAlphasRaw;
    ClAlphas(excludeCl) = [];
    ClRaw = ClData(startCl:stopCl,2);
    Cl = ClRaw;
    Cl(excludeCl) = [];
    
    CdAlphasRaw = CdData(startCd:stopCd,1);
    CdAlphas = CdAlphasRaw;
    CdAlphas(excludeCd) = [];
    CdRaw = CdData(startCd:stopCd,2);
    Cd = CdRaw;
    Cd(excludeCd) = [];
    
    Re = 140000;
end

%% Design Analysis and Plotting

% Plane Constants
SWet = 1.587; % Approximation, m^2
SFuse = 2*0.3*0.08 + 2*0.3*0.2 + 2*pi*(0.04^2) + 2*pi*0.04*0.2 + 2*0.6*0.05 + 2*0.45*0.05;
SRef = 0.397; % Approximation, m^2
Cfe = 0.003; % Civil transport/glider
AR = 14.716; % Aspect Ratio
e = 0.99; % Wing span efficiency

xACWinglet = 0.55; % x position of winglets aerodynamic center

cHoriz = 0.15; % Horizontal stabilizer chord
xACHoriz = 0.825; % x position of horizontal stabilizer aerodynamic center

bVert = 0.2; % Vertical stabilizer chord
xACVert = 0.825; % x position of verticall stabilizer aerodynamic center

xCG = 0.3; % x position of center of gravity
xCGWing = 0.2; % x position of wing CG
xCGFuse = 0.25; % x position of fuselage CG
xCGVert = 0.825; % x position of vertical stab CG
xCGWinglet = 0.55; % x position of the winglet CG
xCGHoriz = 0.825; % x position of horizontal stab CG

xCGPayload = 0.1; % x position of payload CG
WPayload = 0.16; % camera weight, 160 g

xCGBallast = -0.04; % x position of ballast CG
WBallast = 0.05; % Ballast weight

% Other Constants
rhoF = 0.295*9.81; % Foam weight factor, N/m^2
rhoAlt = 1.14; % Air density in Boulder, kg/m^3
h = 17.5; % Start height, m
Vh = 0.6; % Horizontal tail volume coefficient
Vv = 0.035; % Vertical tail volume coefficient
mu = 1.74*(10^-5);

% Analysis

% Find a0
if strcmp(dataSet, dataSet2)
    start = find(ClAlphas <= -7, 1, 'last');
    stop = find(ClAlphas >= 7, 1, 'first');
else
    start = find(ClAlphas <= -5, 1, 'last');
    stop = find(ClAlphas >= 6, 1, 'first');
end
[coef, a0Curve] = leastSquares(ClAlphas(start:stop),Cl(start:stop),1);
a0 = coef(1);

% Find a
a = a0/(1+((57.3*a0)/(pi*e*AR)));

% Find Alpha where L=0

if strcmp(dataSet, dataSet2)
    [~, approxCurve] = leastSquares(ClAlphas,Cl,3);
    alphaL0 = fzero(approxCurve, 0);
else
    [~, approxCurve] = leastSquares(ClAlphas,Cl,5);
    alphaL0 = fzero(approxCurve, -2);
end

% Calculate CL from Cl
CL = a*(ClAlphas - alphaL0);

% Calculate CD from Cd and CL
CD = Cd + ((CL.^2)/(pi*e*AR));

% Calculate full drag polar with Raymer's Oswald factor model
delta = 0.2; % Between 0.2 and 0.3 for low Re

k = 1/(pi*e*AR);
e0 = 0.3/(1 + delta + k*pi*AR);

k1 = 1/(pi*e0*AR)

CDmin = Cfe*(SWet/SRef);

[~, index] = min(CD);

CLminD = CL(index);
k2 = -2*k1*CLminD;

CDo = CDmin + k1*(CLminD)^2

FullCD = CDo + k1*CL.^2 + k2*CL;

Cl_max = max(Cl);
CL_max = max(CL);

%% Performance
% Current
SWingletRaw = 2*(0.5*(0.1+0.15)*(sqrt(0.7088^2 - 0.05^2)));
Sh = 0.03; % 0.3*0.2
SWingletH = 2*(0.5*(0.1+0.15)*(0.4)); % 2 winglets
Sv = 0.03; % 0.3*0.2
SWingletV = 2*(0.5*(0.1+0.15)*(0.2)); % 2 winglets
sweepAngle = atan(2/4);
%SWet = 2*SRef + 2*Sh + 4*Sv + SFuse;
C_L_R = sqrt(CDo./k1)
C_D_R = 2*CDo;
C_L_E = sqrt(3.*CDo./k1)
C_D_E = 4*CDo;
W_span = sqrt(AR.*SRef); % Wing Span
lambda = (1/3).*ones(1,length(W_span));
R_chord = (2.*SRef)./(W_span.*(lambda+1));
MAC = (2/3).*(R_chord).*((1+lambda+lambda.^2)./(1+lambda))

WTO = rhoF*(SRef + SFuse + 2*Sv + Sh + SWingletRaw) + WBallast + WPayload
wingLoadingCurrent = WTO/SRef
LD_Range = C_L_R/C_D_R
LD_Endurance = C_L_E/C_D_E
RmaxCurrent = h*0.5*sqrt(SRef./(k1*Cfe*SWet))
V_RmaxCurrent = sqrt((2.*WTO)./(C_L_R.*rhoAlt.*SRef))
EmaxCurrent = h*(((sqrt((3*CDo)/k1)).^3/2) .* sqrt(0.5*rhoAlt*SRef))./(4*CDo)
V_EmaxCurrent = sqrt((2.*WTO)./(C_L_E.*rhoAlt.*SRef))
xCGCurrent = (rhoF*(0.42*xCGWing + SFuse*xCGFuse + 2*Sv*xCGVert + Sh*xCGHoriz) + WBallast*xCGBallast + WPayload*xCGPayload)./(rhoF*(SRef + SFuse + Sv + Sh) + WBallast + WPayload)
y_mac_trapz = (1/6)*((1+(2*lambda))/(1+lambda))
x_max_trapz = y_mac_trapz*tan(sweepAngle) + 0.25*MAC
y_mac = (W_span/6)*((1+(2*lambda))/(1+lambda))
x_mac = y_mac*tan(sweepAngle) + 0.25*MAC

VhCurrent = (Sh*(xACHoriz-xCGCurrent))/(MAC*SRef) + (SWingletH*(xACWinglet-xCGCurrent))/(MAC*SRef)
VvCurrent = (Sv*(xACVert-xCGCurrent))/(W_span*SRef) + (SWingletV*(xACWinglet-xCGCurrent))/(W_span*SRef)

ARH = (0.2^2)/(0.15*0.2);
a_t = a0/(1+((57.3*a0)/(pi*e*ARH)));
delEpdelAlph = 0.35;
h_np = x_mac/MAC + Vh*(a_t/a)*(1-delEpdelAlph);
hCG = xCGCurrent/MAC;
SM = h_np - hCG

% Varied
SRef = linspace(0,1,100);
W_span = sqrt(AR.*SRef); % Wing Span

xCG = (rhoF*(SRef*xCGWing + SFuse*xCGFuse + Sv*xCGVert + Sh*xCGHoriz) + WBallast*xCGBallast + WPayload*xCGPayload)./(rhoF*(SRef + SFuse + Sv + Sh) + WBallast + WPayload);

Sh = (Vh*SRef*cHoriz)./(xACHoriz - xCG);
Sv = (Vv*SRef.*W_span)./(xACVert - xCG);

SWet = 2*SRef + 2*Sh + 4*Sv + SFuse;

CDo = Cfe*(SWet./SRef);

C_L_R = sqrt(CDo./k1);
C_L_E = sqrt(3.*CDo./k1);

R = h*0.5*sqrt(SRef./(k1*Cfe*SWet));

WTO = rhoF*(SRef + SFuse + 2*Sv + Sh) + WBallast + WPayload;

lambda = (1/3).*ones(1,length(W_span));
R_chord = (2.*SRef)./(W_span.*(lambda+1));
T_chord = lambda.*R_chord;
MAC = (2/3).*(R_chord).*((1+lambda+lambda.^2)./(1+lambda));
CL_Xcr = linspace(0,CL(end),length(W_span));
Re_Cr = 500000;
Xcr = (Re_Cr*mu/rhoAlt).*(sqrt((CL_Xcr.*rhoAlt.*SRef)./(2.*WTO)));

E = h*(((sqrt((3*CDo)/k1)).^3/2) .* sqrt(0.5*rhoAlt*SRef))./(4*CDo);

V_Rmax = sqrt((2.*WTO)./(C_L_R.*rhoAlt.*SRef));
V_Emax = sqrt((2.*WTO)./(C_L_E.*rhoAlt.*SRef));
V_max = sqrt((2.*WTO)./(Cl_max.*rhoAlt.*SRef));

%% Plotting

F = figure();
F.Position = [100 100 940 740];

sgtitle("Drag Polar and Lift Curve for Design")

% Cl/CL vs. alpha
subplot(1,3,1)
hold on;
grid on;

AlphaCl2D = plot(ClAlphasRaw, ClRaw);
AlphaCL = plot(ClAlphas, CL);

% Utility lines
% alphaTest = -5:0.001:5;
% plot(alphaTest, approxCurve(alphaTest));
% plot(ClAlphas, a0Curve(ClAlphas));
alpha0Line = xline(alphaL0,'m--');
alpha0Label = sprintf("\\alpha_{L=0} = %.3f^o", alphaL0);
xline(0);
yline(0);

% Title, legend, labels
subset = [AlphaCl2D, AlphaCL, alpha0Line];
titles = ["\alpha vs. C_l", "\alpha vs. C_L, calculated", alpha0Label];

title('Lift Coefficient vs. Angle of Attack')
xlabel('\alpha (deg)')
ylabel('Coefficient of lift')
legend(subset, titles, 'Location', 'best');
hold off;


% Cd/CD vs. alpha
subplot(1,3,2)
if strcmp(dataSet, dataSet2)
    hold on;
    grid on;
    
    AlphaCd2D = plot(CdAlphasRaw, CdRaw);
    AlphaCD = plot(CdAlphas, CD);
    
    % Utility lines
    % alphaTest = -5:0.001:5;
    % plot(alphaTest, approxCurve(alphaTest));
    % plot(ClAlphas, a0Curve(ClAlphas));
    % alpha0Line = xline(alphaL0,'m--');
    % alpha0Label = sprintf("\\alpha_{L=0} = %.3f^o", alphaL0);
    xline(0);
    yline(0);
    
    % Title, legend, labels
    subset = [AlphaCd2D, AlphaCD];
    titles = ["\alpha vs. C_d", "\alpha vs. C_D, calculated"];
    
    title('Drag Coefficient vs. Angle of Attack')
    xlabel('\alpha (deg)')
    ylabel('Coefficient of drag')
    legend(subset, titles, 'Location', 'best');
    hold off;
end

% Drag polar
subplot(1,3,3)
grid on;
hold on;

% startCl = find(ClAlphas <= -10, 1, 'last')
% stopCl = find(ClAlphas >= 10, 1, 'first')
% 
% startCd = find(CdAlphas <= -10, 1, 'last')
% stopCd = find(CdAlphas >= 10, 1, 'first')

DragPolar2D = plot(Cl, Cd);
DragPolar3D = plot(CL, CD);
FullDragPolar = plot(CL, FullCD);

% Utility lines
xline(0);
yline(0);

% Title, legend, labels
subset = [DragPolar2D, DragPolar3D, FullDragPolar];
titles = ["C_d vs. C_l", "C_D vs. C_L, calculated", "Full Aircraft Drag Polar"];

title('Drag Polar')
xlabel('Coefficient of lift')
ylabel('Coefficient of drag')
legend(subset, titles, 'Location', 'best')

hold off;


G = figure();
G.Position = [100 100 940 740];

sgtitle("Design performance")

subplot(2,3,1)
hold on;
title('Max R & Cl Req vs W/S')
xlabel('Wing Loading Required - W/Sref(N/m^2)')
plot(WTO./SRef, R)
yline(125,'--b')
yline(175,'-.b')
yyaxis left
ylabel('Range (m)')
yyaxis right
plot(WTO./SRef, C_L_R)
yline(CL_max,'--r')
ylabel('Coefficient of Lift Required')
legend('R Achieved','Rmin Req','Rmax Req','CL Req for R','CLmax Limit')
xlim([10 80])
hold off

subplot(2,3,2)
hold on;
plot(WTO./SRef, E)
yline(13,'--b')
yline(20,'-.b')
title('Max E & Cl Req vs W/S')
xlabel('Wing Loading Required - W/Sref(N/m^2)')
yyaxis left
ylabel('Endurance Achieved (s)')
yyaxis right
plot(WTO./SRef, C_L_E)
yline(CL_max,'--r')
ylabel('Coefficient of Lift Required')
legend('E Achieved','Emin Req','Emax Req','CL Req for E','CLmax Limit')
xlim([10 80])
hold off

subplot(2,3,3)
hold on;
plot(WTO./SRef, WTO);
title('Variation of Aircraft Weight with Wing Loading')
xlabel('Wing Loading - W/Sref(N/m^2)')
ylabel('Aircraft Weight - Wto (N)')
xlim([10 80])

subplot(2,3,4)
hold on
title('Vel Req for Max R & Max E')
xlabel('Wing Loading Required - W/Sref(N/m^2)')
yline(7,'--g')
yline(12,'-.g')
plot(WTO./SRef, V_Emax)
plot(WTO./SRef, V_Rmax,'g')
plot(WTO./SRef, V_max,'--r')
ylabel('Velocity Required (m/s)')
xlim([10 80])
legend('Vmin for R','Vmax for R','Vel to Achieve Emax','Vel to Achieve Rmax','Limit Vel for CLmax')
hold off

subplot(2,3,5)
hold on
title('Variation of Geometry with Constant AR, Taper, VH, VV, lt, lv')
xlabel('Wing Loading - W/Sref(N/m^2)')
yyaxis left
plot(WTO./SRef,W_span)
plot(WTO./SRef,MAC,'--')
plot(WTO./SRef,R_chord,':')
plot(WTO./SRef,T_chord,'-.')
plot(WTO./SRef,Xcr,'--k')
ylabel('Length(m)')
yyaxis right
plot(WTO./SRef,SRef)
plot(WTO./SRef,Sh,'--')
plot(WTO./SRef,Sv,':')
ylabel('Planform Area - Sref (m^2)')
xlim([10 80])
legend('Wing Span','MAC','Root Chord','Tip Chord','Laminar Trans Xcr','Sref wing','Sref horizontal stabilizer','Sref vertical stabilizer')
hold off

subplot(2,3,6)
hold on;
plot(WTO./SRef, 100*(WTO-WPayload))
title('Variation of Aircraft Cost with Wing Loading')
xlabel('Wing Loading - W/Sref(N/m^2)')
ylabel('Total Aircraft Cost ($)')
xlim([10 80])

%-------------------------------------------------------------------------%

%% Least Squares Function from ASEN 2012

function [X,f] = leastSquares(t,y,p)
    % for writing this function, some skeleton code has been provided to
    % help you design the function to serve your purposes
    A = [];
    % write an expression for A, the input matrix
    for ii = 0:p
        col = t.^ii;
        A = [col, A];
    end
    % compute coefficient vector, x_hat
    x_hat = A\y;
    X = x_hat;
    
    % do not change the following lines of code. This will generate the
    % anonymous function handle "f" for you
%     f = '@(x)';
%     for i = 0:p
%         f = strcat(f,'+',strcat(string(x_hat(i+1)),'.*x.^',string(p-i)));
%     end
%     eval(strcat('f = ',f,';'))
    
    while length(x_hat) < 7
        x_hat = [0;x_hat];
    end
    % workaround for MATLAB grader
    f = @(x) x_hat(1)*x.^6 + x_hat(2)*x.^5 + x_hat(3)*x.^4 + x_hat(4)*x.^3 + x_hat(5)*x.^2 + x_hat(6)*x + x_hat(7);
    
end






