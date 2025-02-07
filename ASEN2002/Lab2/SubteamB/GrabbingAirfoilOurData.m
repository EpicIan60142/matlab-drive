clear;
clc;

numFiles = 6;

%% file loading our files
for i= 1:numFiles %find name of each group data file
    % grab the data name and put it in a string
    files{i} = dir(['AirfoilPressure_S2_',num2str(i),'*.csv']); 
    % maybe name for each section like "files1030" 
end


for j=1:numFiles % load in each test data 
  % load data using the data name into an array
  test = load(files{j}.name);    % need to call info{1,#} to get data out
  

 % find index where angle of attack change and use them to make smaller sections
 % of data
 % create a vector of angle of attack
 AOA = test(:,7);
 
 %create a vector of velocities
 vel = test(:,4); 
 
 % create a vector of indexes where aoa changes or velocity changes
 AOA_Change_Index = find(abs(diff(AOA)) > 0.5);
 vel_change_Index = find(diff(vel) > 2);
 
 % add start and end to the indexes -- these are the different groups of
 % aoa
 AOAVel_Change_Index = [0;AOA_Change_Index;vel_change_Index;length(AOA)];
 AOAVel_Change_Index = sort(AOAVel_Change_Index);
 
 % create a vector of indexes where velocity changes 
 %  vel_Change_Index = find(diff(vel) > 2);
 %  vel_Change_Index = [0;vel_Change_Index;length(vel)];
 
 % preallocate array for speed
 sections = zeros(length(AOAVel_Change_Index)-1,29);
 
 
 % loop through and put data into multiple sections
 for k = 1:length(AOAVel_Change_Index)-1
     % collapse the data into a row for each section by taking mean of the
     % section
     sections(k,:)= mean(test((AOAVel_Change_Index(k)+1:AOAVel_Change_Index(k+1)),:));
 end 
 
 % take each angle of attack and further split it into sections based on
 % velocity
 
 
%% need to put measurements in structs  

    % putting all the measurement groups in vectors 
    atmTemp = sections(:,1);
    atmPress = sections(:,2);
    atmRho = sections(:,3);
    freestreamV = sections(:,4);
    pitotDynPress1 = sections(:,5);
    pitotDynPress2 = sections(:,6);
    AngAttack = sections(:,7);
    Normal = sections(:,8);
    Axial = sections(:,9);
    Pitch = sections(:,10);
    PressPort = [];
    
    for p = 1:16
    PressPort(:,p) = sections(:,p+13);
    end

    % put each group into structs 
    AirfoilInfo(j).atmPress = atmPress;
    AirfoilInfo(j).atmTemp = atmTemp;
    AirfoilInfo(j).atmRho = atmRho; 
    AirfoilInfo(j).freestreamV = freestreamV;
    AirfoilInfo(j).pitotDynPress1 = pitotDynPress1;
    AirfoilInfo(j).pitotDynPress2 = pitotDynPress2;
    for p = 1:16
        AirfoilInfo(j).PressPort(:,p) = PressPort(:,p);
    end
    AirfoilInfo(j).AngAttack = AngAttack;
    AirfoilInfo(j).Normal = Normal;
    AirfoilInfo(j).Axial = Axial;
    AirfoilInfo(j).Pitch = Pitch;    

end

AllData = struct2table(AirfoilInfo);

% Combining all the 

Final.atmPress = [];
Final.atmTemp = [];
Final.atmRho = [];
Final.freestreamV = [];
Final.pitotDynPress1 = [];
Final.pitotDynPress2 = [];
Final.PressPort1 = [];
Final.PressPort2 = [];
Final.PressPort3 = [];
Final.PressPort4 = [];
Final.PressPort5 = [];
Final.PressPort6 = [];
Final.PressPort7 = [];
Final.PressPort8 = [];
Final.PressPort9 = [];
Final.PressPort10 = [];
Final.PressPort11 = [];
Final.PressPort12 = [];
Final.PressPort13 = [];
Final.PressPort14 = [];
Final.PressPort15 = [];
Final.PressPort16 = [];
Final.AngAttack = [];
Final.Normal = [];
Final.Axial = [];
Final.Pitch = [];
PPtable1 = [];
PPtable2 =[];
PPtable3 = [];
PPtable4 = [];
PPtable5 = [];
PPtable6 = [];
PPtable7 = [];
PPtable8 = [];
PPtable9 = [];
PPtable10 = [];
PPtable11 = [];
PPtable12 = [];
PPtable13 = [];
PPtable14 = [];
PPtable15 = [];
PPtable16 = [];

for k=1:numFiles %length(AllData.atmPress) % number of num files
Final.atmPress = [Final.atmPress; AllData.atmPress{k}];
Final.atmTemp = [Final.atmTemp; AllData.atmTemp{k}];
Final.atmRho = [Final.atmRho; AllData.atmRho{k}];
Final.freestreamV = [Final.freestreamV; AllData.freestreamV{k}];
Final.pitotDynPress1 = [Final.pitotDynPress1; AllData.pitotDynPress1{k}];
Final.pitotDynPress2 = [Final.pitotDynPress2; AllData.pitotDynPress2{k}];
  
PPtable1 = [PPtable1; AllData.PressPort{k}(:,1)];
PPtable2 = [PPtable2; AllData.PressPort{k}(:,2)];
PPtable3 = [PPtable3; AllData.PressPort{k}(:,3)];
PPtable4 = [PPtable4; AllData.PressPort{k}(:,4)];
PPtable5 = [PPtable5; AllData.PressPort{k}(:,5)];
PPtable6 = [PPtable6; AllData.PressPort{k}(:,6)];
PPtable7 = [PPtable7; AllData.PressPort{k}(:,7)];
PPtable8 = [PPtable8; AllData.PressPort{k}(:,8)];
PPtable9 = [PPtable9; AllData.PressPort{k}(:,9)];
PPtable10 = [PPtable11; AllData.PressPort{k}(:,10)];
PPtable11 = [PPtable11; AllData.PressPort{k}(:,11)];
PPtable12 = [PPtable12; AllData.PressPort{k}(:,12)];
PPtable13 = [PPtable13; AllData.PressPort{k}(:,13)];
PPtable14 = [PPtable14; AllData.PressPort{k}(:,14)];
PPtable15 = [PPtable15; AllData.PressPort{k}(:,15)];
PPtable16 = [PPtable16; AllData.PressPort{k}(:,16)];

Final.PressPort1 = PPtable1;
Final.PressPort2 = PPtable2;
Final.PressPort3 = PPtable3;
Final.PressPort4 = PPtable4;
Final.PressPort5 = PPtable5;
Final.PressPort6 = PPtable6;
Final.PressPort7 = PPtable7;
Final.PressPort8 = PPtable8;
Final.PressPort9 = PPtable9;
Final.PressPort10 = PPtable10;
Final.PressPort11 = PPtable11;
Final.PressPort12 = PPtable12;
Final.PressPort13 = PPtable13;
Final.PressPort14 = PPtable14;
Final.PressPort15 = PPtable15;
Final.PressPort16 = PPtable16;

Final.AngAttack = [Final.AngAttack; AllData.AngAttack{k}];
Final.Normal = [Final.Normal; AllData.Normal{k}];
Final.Axial = [Final.Axial; AllData.Axial{k}];
Final.Pitch = [Final.Pitch; AllData.Pitch{k}];
end

%% Find Airspeed calculations 
const.chordL = 3.5/39.37; % in to m



