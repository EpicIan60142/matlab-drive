clear
clc

ind = 0;

h_drop = 80.5358;

for i = 0:170
   h =  h_drop * (1/(1+exp((i-70)/8)));
   
   ind = ind +  1;
   
   xvals(ind) = ind;
   
   hvals(ind) = h;
end

figure
plot(xvals,hvals)

