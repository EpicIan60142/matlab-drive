function slope = eulerFunc(t,A,p)
%eulerFunc: Gives the slope of a differential equation for further
%approximation
%Inputs: t = time (years), A = amount owed, p = monthly payment

if(t <= 5)
    slope = 0.03*A - 12*p;
elseif(t > 5)
    slope = (0.03+0.015*sqrt(t-5))*A - 12*p;
end

end

