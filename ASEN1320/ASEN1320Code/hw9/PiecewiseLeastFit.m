function Yvector = PiecewiseLeastFit(M,B,Xvector)

N = length(Xvector);
Yvector = zeros(1,N);

for n = 1 : N
   
    x = Xvector(n);
   
    if x < 10
       Yvector(n) = M(1) * Xvector(n) + B(1);
    elseif x >= 10  && x < 15
       Yvector(n) = M(2) * Xvector(n) + B(2);
    elseif x >= 15  && x < 20
       Yvector(n) = M(3) * Xvector(n) + B(3);
    else % x >= 20  
       Yvector(n) = M(4) * Xvector(n) + B(4);
    end

end

end