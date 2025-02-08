function input = extractData(fileName, rho, c, k, alpha, d)
% extractData - Pulls in a file name and various physical constants to
%               create the base input structure for all ASEN 3113 Lab 2 
%               functions

    input = struct('matTitle',[],'matProp',[],'time',[],'thermocouples',[], ...
                    'x',[],'qDot',[],'area',[]);

    data = readmatrix(fileName);
    
    input.time = data(:,1);
    if size(data,2) == 10 % T0 included
        T0 = data(:,2);
        T1 = data(:,3);
        T2 = data(:,4);
        T3 = data(:,5);
        T4 = data(:,6);
        T5 = data(:,7);
        T6 = data(:,8);
        T7 = data(:,9);
        T8 = data(:,10);
        input.thermocouples = [T0, T1, T2, T3, T4, T5, T6, T7, T8];
        input.x = [0, 1.375:0.5:1.375+7*0.5];
    else % No T0
        T1 = data(:,2);
        T2 = data(:,3);
        T3 = data(:,4);
        T4 = data(:,5);
        T5 = data(:,6);
        T6 = data(:,7);
        T7 = data(:,8);
        T8 = data(:,9);
        input.thermocouples = [T1, T2, T3, T4, T5, T6, T7, T8];
        input.x = 1.375:0.5:1.375+7*0.5;
    end
    
    content = extractAfter(fileName,"Data\");
    material = extractBefore(content,"_");
    volt = str2double(erase(extractBetween(content,"_","_"),'V')); % V
    current = str2double(erase(extractAfter(content,"V_"),'mA'))/1000; % A
    
    input.matTitle = sprintf("%s at %.0f V and %.0f mA", material, volt, current*1000);

    switch material
        case "Aluminum"
            input.matProp = [rho(1), c(1), k(1), alpha(1)];
        case "Brass"
            input.matProp = [rho(2), c(2), k(2), alpha(2)];
        case "Steel"
            input.matProp = [rho(3), c(3), k(3), alpha(3)];
        otherwise
    end
    
    input.qDot = volt*current; % J
    
    input.area = pi*(d/2)^2; % m^2

end

