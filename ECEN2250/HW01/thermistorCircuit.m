%% ECEN2250 HW01, Ian Faber, 08/29/2021
clc;clear;close all;

%% Problem 3c: Thermistor Sim
Vdd = 3.3; %Voltage regulator voltage
NTCResistance = [5774,7097]; %Range of NTC resistances
resistance = 30000; %resistance of an SMD resistor
costSingle = 0.10; %Cost of one SMD resistor ($0.10)
costTen = 0.12; %Cost of ten SMD resistors ($0.12)
maxCurrent = 10e-6; %Maximum current allowed for the circuit in Amps

fprintf("Simulation with minimum NTC Resistance (5774 Ohms):\n\n")
for k = 1:11 %Loop from 1 30000 ohm resistor to 10 in series
    if(k==10)
        cost = costTen;
    elseif(k > 10)
        cost = costTen + (k-10)*costSingle;
    else
        cost = k*costSingle;
    end
    resistanceT = NTCResistance(1) + k*resistance; 
    current = Vdd / resistanceT; %Ohm's Law
    if(current > maxCurrent)
        fprintf("Circuit NOT COMPLIANT with %d resistor(s), max current is %f (> 10) uA, but the cost would be $%f\n",k,current*1e6,cost);
    else
        fprintf("Circuit IS COMPLIANT with %d resistor(s), max current is %f uA and cost is $%f\n",k,current*1e6,cost);
    end
end
fprintf("\nSimulation with maximum NTC Resistance (7097 Ohms):\n\n")
for k = 1:11 %Loop from 1 30000 ohm resistor to 10 in series
    if(k==10)
        cost = costTen;
    elseif(k > 10)
        cost = costTen + (k-10)*costSingle;
    else
        cost = k*costSingle;
    end
    resistanceT = NTCResistance(2) + k*resistance; 
    current = Vdd / resistanceT; %Ohm's Law
    if(current > maxCurrent)
        fprintf("Circuit NOT COMPLIANT with %d resistor(s), max current is %f (> 10) uA, but the cost would be $%f\n",k,current*1e6,cost);
    else
        fprintf("Circuit IS COMPLIANT with %d resistor(s), max current is %f uA and cost is $%f\n",k,current*1e6,cost);
    end
end
    