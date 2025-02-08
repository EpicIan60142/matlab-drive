function dfdt = finiteDifMat(t0, dt, f)
% Caculates a numerical difference quotient for a 3x3 matrix

f1 = cell2mat(f(cell2mat(f(:,1)) == t0+dt, 2));
f2 = cell2mat(f(cell2mat(f(:,1)) == t0, 2));
dfdt = (f1 - f2)./dt;

end