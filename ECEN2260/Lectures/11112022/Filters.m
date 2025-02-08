% Housekeeping
clc; clear; close all

% Cut-off Frequency 
fc = 1000;      % Hz
wc = 2*pi*fc;   % rad/s

% Ripple Band (for Chebychev)
rb = 3;         % dB

%% Butterworth Filter
for n = 1:5
    % Compute coefficients of the numerator (b) and denominator (a)
    [b,a] = butter(n,wc,'s');
    
    % Define the transfer function
    G = tf(b,a);
        
    % Show pole locations on the complex plane
    figure(1)
    pzmap(G)
	hold on
    grid on
    
    % Draw the Bode Diagram
    figure(2)
    bode(G,{wc*1e-2,wc*1e2})
    hold on
    grid on

end

figure(1)
legend('n=1','n=2','n=3','n=4','n=5')
title('Pole-Zero Map: Butterworth Filter')

figure(2)
legend('n=1','n=2','n=3','n=4','n=5')
title('Bode Diagram: Butterworth Filter')

%% Chebychev Filter
for n = 1:5
    % Compute coefficients of the numerator (b) and denominator (a)
    [b,a] = cheby1(n,rb,wc,'s');
    
    % Define the transfer function
    H = tf(b,a);
        
    % Show pole locations on the complex plane
    figure(3)
    pzmap(H)
	hold on
    grid on
    
    % Draw the Bode Diagram
    figure(4)
    bode(H,{wc*1e-2,wc*1e2})
    hold on
    grid on
end

figure(3)
legend('n=1','n=2','n=3','n=4','n=5')
title('Pole-Zero Map: Chebychev Filter')

figure(4)
legend('n=1','n=2','n=3','n=4','n=5')
title('Bode Diagram: Chebychev Filter')

%% Comparison between the two filters for n=5
figure
bode(G,{wc*1e-2,wc*1e2})
hold on
bode(H,{wc*1e-2,wc*1e2})
grid on
legend('Butterworth Filter','Chebychev Filter')
title('Bode Diagram: Filter comparison for n=5')
