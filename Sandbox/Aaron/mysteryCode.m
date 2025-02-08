clc;
clear all;
close all;
%% syracuse problem try 1

%% create radom table of starting values
nt   = 100; % number of trials
rtmi = 1; % min numbers of the random numbers
rtma = 100; % max numbers of the random numbers

startm = 1:nt;
for i = 1:nt
    stnum = randi([rtmi rtma]);
    startm(i) = stnum;
end

bigm = ones(nt);     % creating larger matrix for preallocation of values
bigm(1,:) = startm'; % allocating the first random numbers generated into the larger matrix

%% Make loop that both finds new values from orginal values and stores the values in a matrix
 
for j = 1:nt
    count = 1;
  while bigm(count,j)> 1 % checks if its over one since loop repeats after this
    if   bigm(count,j) == 1 % checks if its over one  
         bigm(count:end,j) = 1;% and if so sets the rest of the values equal to 0 as well
    elseif mod(bigm(count,j),2) == 0 % checks if the number previolsy found has a reaminder equal to 0
       % since this would imply the numner is even
        bigm(count + 1,j) =  bigm(count,j)/2;  
    else
        % should be odd numbers since code has already elimanted all even
        % ones
        bigm(count + 1,j) = bigm(count,j)*3 +1 ;
    end
     count = count + 1; % keeps count to itertate
  end
 
  T = length(bigm); % finds larger dimension of the matrix which repersents the needed max value on the x axis
 
  % ploting the count vs the number at each count
  figure(1)
  plot((1:T),bigm(:,j))
  hold on
  title('Syracuse Problem (count vs number)','FontSize',16);
  xlabel('Count, m','FontSize',16);
  ylabel('Numerical value, m','FontSize',16);
 
end

%% Start of anayalsis

sam  = zeros(2,nt); % preallocating matrix for anyalsyis

% upcoming for loop helps find how long it took each starting number to get to  value '1'
for k = 1:nt
    sam(1,k) = bigm(1,k);
    Ta1 = find(bigm(:,k) > 1);
    sam(2,k) = length(Ta1);
end

% upcoming graph helps to display the retaionship between the starting number and how long it takes to get to '1'
    figure(2)
    scatter(sam(2,:),sam(1,:))
    hold on
    title('Syracuse Problem (count vs starting number)','FontSize',16);
    xlabel('Count, m','FontSize',16);
    ylabel('Starting numerical value, m','FontSize',16);
