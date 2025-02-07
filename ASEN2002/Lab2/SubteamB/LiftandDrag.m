clear;
clc;
close all;
% loads in our final strtuct and port locations
load('Data.mat');
ClarkY14NACA = readtable('ClarkY14_NACA_TR628.xlsx','Range','A2:C12');

% put how many different combinations of angles of attack paired with
% different velcoties there are
states = length(Final.AngAttack);
Cps = zeros(states,17);

% loops through each state
for State = 1:states
Cp = ones(1,16);

% grab the freestream variables from the final struct
Pinf = Final.atmPress(State);
Rhoinf = Final.atmRho(State);
Vinf = Final.freestreamV(State);

% grab each pressure port pressure    
    P(1) = Final.PressPort1(State);
    P(2) = Final.PressPort2(State);
    P(3) = Final.PressPort3(State);
    P(4) = Final.PressPort4(State);
    P(5) = Final.PressPort5(State);
    P(6) = Final.PressPort6(State);
    P(7) = Final.PressPort7(State);
    P(8) = Final.PressPort8(State);
    P(9) = Final.PressPort9(State);
    P(10) = Final.PressPort10(State);
    P(11) = Final.PressPort11(State);
    P(12) = Final.PressPort12(State);
    P(13) = Final.PressPort13(State);
    P(14) = Final.PressPort14(State);
    P(15) = Final.PressPort15(State);
    P(16) = Final.PressPort16(State);
    
    Cp_no_approx = ones(1,16);
    for i=1:16
        % calculate a coefficeint for each pressure port
        Cp_no_approx(i) = (P(i))/(0.5*Rhoinf*Vinf^2);
    end
    % estimate the trailing edge coefficient of pressure by using our
    % function
    Cp_trailing_edge = Estimate_Trailing_Edge(Cp_no_approx,Port_Locations);
    % add in the trailing edge coefficient in it's placement relative to
    % the other pressure ports
    Cp_approx = [Cp_no_approx(1:9),Cp_trailing_edge,Cp_no_approx(10:end)];
    
    % for each State - put the coeffiecients of pressure in a row
    Cps(State,:) = Cp_approx;
end

% find the chord
chord = max(Port_Locations(1,:));
% normalize the locations relative to the chord -- easier to do
% calculations
normalized_Locations = Port_Locations/chord;

% plot the normalized locations of the pressure ports versus the coeffients
% of pressure for angle of attack 10 and vinf 33m/s -- the last State



% calling Integration to get lift and drag coefficients
[Cl, Cd] = Integration(Cps, normalized_Locations, Final);

% need to split the angles of attack into the different velocites since
% there are 3 different velocities per angle of attack
% We must analyze them seperately or it will skew results 
counter = 1;
for i = 1:3:(length(Final.AngAttack))
    Cl_10(counter) = Cl(i);
    Cl_18(counter) = Cl(i+1);
    Cl_33(counter) = Cl(i+2);
    Cd_10(counter) = Cd(i);
    Cd_18(counter) = Cd(i+1);
    Cd_33(counter) = Cd(i+2);
    Cps_10(counter,:) = Cps(i,:);
    Cps_18(counter,:) = Cps(i+1,:);
    Cps_33(counter,:) = Cps(i+2,:);

    counter = counter+1;
end

AoA = Final.AngAttack(1:3:end)';

%% Sort the stuff
%% Sort Cp
%for velocity 10m/s
Cps_10_AoA = [AoA', Cps_10];
% sort in ascending order according to angle of attack
[~,index] = sort(Cps_10_AoA(:,1));
sort_Cps_10_AoA = Cps_10_AoA(index,:);

%for velocity 18m/s
Cps_18_AoA = [AoA', Cps_18];
% sort in ascending order according to angle of attack
[~,index] = sort(Cps_18_AoA(:,1));
sort_Cps_18_AoA = Cps_18_AoA(index,:);

%for velocity 33m/s
Cps_33_AoA = [AoA', Cps_33];
% sort in ascending order according to angle of attack
[~,index] = sort(Cps_33_AoA(:,1));
sort_Cps_33_AoA = Cps_33_AoA(index,:);
%% Sort Cl 

%for velocity 10m/s
Cl_10_AoA = [AoA; Cl_10]';
% sort in ascending order according to angle of attack
[~,index] = sort(Cl_10_AoA(:,1));
sort_Cl_10_AoA = Cl_10_AoA(index,:);

%for velocity 18m/s
Cl_18_AoA = [AoA; Cl_18]';
% sort in ascending order according to angle of attack
[~,index] = sort(Cl_18_AoA(:,1));
sort_Cl_18_AoA = Cl_18_AoA(index,:);

%for velocity 33m/s
Cl_33_AoA = [AoA; Cl_33]';
% sort in ascending order according to angle of attack
[~,index] = sort(Cl_33_AoA(:,1));
sort_Cl_33_AoA = Cl_33_AoA(index,:);

%% Sort Cd

%for velocity 10m/s
Cd_10_AoA = [AoA; Cd_10]';
% sort in ascending order according to angle of attack
[~,index] = sort(Cd_10_AoA(:,1));
sort_Cd_10_AoA = Cd_10_AoA(index,:);

%for velocity 18m/s
Cd_18_AoA = [AoA; Cd_18]';
% sort in ascending order according to angle of attack
[~,index] = sort(Cd_18_AoA(:,1));
sort_Cd_18_AoA = Cd_18_AoA(index,:);

%for velocity 33m/s
Cd_33_AoA = [AoA; Cd_33]';
% sort in ascending order according to angle of attack
[~,index] = sort(Cd_33_AoA(:,1));
sort_Cd_33_AoA = Cd_33_AoA(index,:);
 

%% plotting
% plot angle of attack vs coeffient of lift with velcity 10 
figure(1); hold on;
%Cps for different AoA's
%33 m/s
%AoA = -9, 2, 7, 15
Cps_we_want_to_plot_fig1 = [sort_Cps_33_AoA(5,2:end);sort_Cps_33_AoA(15,2:end);sort_Cps_33_AoA(18,2:end);sort_Cps_33_AoA(24,2:end)];
grid on; grid minor;
plot(normalized_Locations(1,:), Cps_we_want_to_plot_fig1,'.-','MarkerSize',15)
set(gca, 'YDir','reverse')
xlabel('Normalized Locations (chord)')
ylabel('Coefficients of Pressure')
title({'Coefficients of Pressure vs Normalized Locations','for Selected Angles of Attack at 33m/s'})
leg = legend('-9','2','7','15','location','northeast');
title(leg,'Angles of Attack');
saveas(figure(1),'Cps vs AoA.png')
hold off;


figure(2); hold on;
grid on; grid minor;
% plot(sort_Cl_10_AoA(:,1), 0.80*sort_Cl_10_AoA(:,2),'.-r','MarkerSize',15);
% plot(sort_Cl_18_AoA(:,1), 0.80*sort_Cl_18_AoA(:,2),'.-b','MarkerSize',15);
% plot(sort_Cl_33_AoA(:,1), 0.80*sort_Cl_33_AoA(:,2),'.-g','MarkerSize',15);
plot(sort_Cl_10_AoA(:,1), sort_Cl_10_AoA(:,2),'.-r','MarkerSize',15);
plot(sort_Cl_18_AoA(:,1), sort_Cl_18_AoA(:,2),'.-b','MarkerSize',15);
plot(sort_Cl_33_AoA(:,1), sort_Cl_33_AoA(:,2),'.-g','MarkerSize',15);
plot(table2array(ClarkY14NACA(:,1)), table2array(ClarkY14NACA(:,2)),'.-k','MarkerSize',15);
xlabel('Angle of Attack (°)')
ylabel('Coefficents of Lift')
title('Coefficients of Lift vs Angle of Attack for Different Velocities')
legend('Velocity 10 m/s','Velocity 18 m/s','Velocity 33 m/s','Clark Y14 NACA Data','Location','northwest')
saveas(figure(2),'Cls vs AoA.png')
hold off;

% plot angle of attack vs coeffient of lift with velcity 10 
figure(3); hold on;
grid on; grid minor;
plot(sort_Cd_10_AoA(:,1), sort_Cd_10_AoA(:,2),'.-r','MarkerSize',15);
plot(sort_Cd_18_AoA(:,1), sort_Cd_18_AoA(:,2),'.-b','MarkerSize',15);
plot(sort_Cd_33_AoA(:,1), sort_Cd_33_AoA(:,2),'.-g','MarkerSize',15);
plot(table2array(ClarkY14NACA(:,1)), table2array(ClarkY14NACA(:,3)),'.-k','MarkerSize',15);
xlabel('Angle of Attack (°)')
ylabel('Coefficents of Drag')
title('Coefficients of Drag vs Angle of Attack for Different Velocities')
legend('Velocity 10 m/s','Velocity 18 m/s','Velocity 33 m/s','Clark Y14 NACA Data','Location','northwest')
saveas(figure(3),'Cds vs AoA.png')
hold off;

% plotting cd and cl on same graph
figure(4); hold on;
grid on; grid minor;

yyaxis left
plot(sort_Cl_10_AoA(:,1), sort_Cl_10_AoA(:,2),'o-r','MarkerSize',10);
plot(sort_Cl_18_AoA(:,1), sort_Cl_18_AoA(:,2),'s-r','MarkerSize',10);
plot(sort_Cl_33_AoA(:,1), sort_Cl_33_AoA(:,2),'^-r','MarkerSize',10);
plot(table2array(ClarkY14NACA(:,1)), table2array(ClarkY14NACA(:,2)),'.-r','MarkerSize',15);
xlabel('Angle of Attack (°)')
ylabel('Coefficents of Lift')
% title('Angle of Attack vs Coefficients of Lift for Different Velocities')
% legend('Velocity 10 m/s','Velocity 18 m/s','Velocity 33 m/s','Location','northwest')

yyaxis right
plot(sort_Cd_10_AoA(:,1), sort_Cd_10_AoA(:,2),'o-b','MarkerSize',7);
plot(sort_Cd_18_AoA(:,1), sort_Cd_18_AoA(:,2),'s-b','MarkerSize',7);
plot(sort_Cd_33_AoA(:,1), sort_Cd_33_AoA(:,2),'^-b','MarkerSize',7);
plot(table2array(ClarkY14NACA(:,1)), table2array(ClarkY14NACA(:,3)),'.-b','MarkerSize',15);

ylabel('Coefficents of Drag')

title('Coefficients of Lift and Drag vs Angle of Attack for Different Velocities')
legend('Cl Velocity 10 m/s','Cl Velocity 18 m/s','Cl Velocity 33 m/s',...
    'Cl Clark Y14 NACA','Cd Velocity 10 m/s','Cd Velocity 18 m/s','Cd Velocity 33 m/s',...
    'Cd Clark Y14 NACA','Location','northwest')
saveas(figure(4),'Cls and Cds vs AoA.png')
hold off;

figure(5); hold on;
Cps_we_want_to_plot = [sort_Cps_10_AoA(18,2:18);sort_Cps_18_AoA(18,2:18);sort_Cps_33_AoA(18,2:18)];
grid on; grid minor;
plot(normalized_Locations(1,:), Cps_we_want_to_plot,'.-','MarkerSize',15)
set(gca, 'YDir','reverse')
xlabel('Normalized Locations (chord)')
ylabel('Coefficients of Pressure')
title({'Coefficients of Pressure vs Normalized', 'Locations for 7° Angle of Attack'})
leg = legend('10 m/s','18 m/s','33 m/s','location','southeast');
title(leg,'Velocities');
saveas(figure(5),'Cps vs velocity.png')
hold off;

