% Write a scrip (mainFit) and three functions (GenerateData, LeastSquares, PiecewiseLeastFit)

xData = linspace(-5,25); 
yData = GenerateData(xData); % Check the output from GenerateData (3pt)

M = zeros(1,4); 
B = zeros(1,4); 
ind = xData < 10;
[M(1),B(1)] = LeastSquares(xData(ind),yData(ind)); %Check the output from LeastSquares (5pt)
ind = (xData < 15 & xData >= 10);
[M(2),B(2)] = LeastSquares(xData(ind),yData(ind)); 
ind = xData < 20 & xData >= 15;
[M(3),B(3)] = LeastSquares(xData(ind),yData(ind)); 
ind = xData >= 20;
[M(4),B(4)] = LeastSquares(xData(ind),yData(ind)); 

%Is it possible to check on M and B (3pt)?

xFit = linspace(-5,25,1000);
yFit  = PiecewiseLeastFit(M,B,xFit); %Check this output (5pt)

%Ask students to plot the output for their testing but no need to check
plot(xData, yData,' o'); 
hold on
plot(xFit, yFit);
xlabel('Angle of attack, \alpha (in degrees)');
ylabel('Coefficient of lift, C_l');
legend('Experimental Data', 'Piecewise Linear Fit');
ax=gca; ax.FontSize= 16;ax.LineWidth = 1;
