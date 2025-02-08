function pumpkinPlot()
%HALLOWEENGRAPH makes a halloween graph
%   no output
%   input: days array the day of Halloween month
%          shape 1 = ghost, 2 = bat, 3 = pumpkin

r = pi*(-24:1:24)/24; 
s = pi*(-24:1:24)/24; 
[theta phi] = meshgrid(r,s); 
f = '(2 + cos(phi))*cos(theta)'; 
g = '(2 + cos(phi))*sin(theta)'; 
h = 'sin(phi)'; 
F = vectorize(f); 
G = vectorize(g); 
H = vectorize(h); 
x = eval(F); 
y = eval(G); 
z = eval(H); 
hSurface = surf(x,y,z);
set(hSurface,'FaceColor',[0.8500 0.3250 0.0980])
hold on
theta2 = (0:1:36)*pi/18;
rho2 = [0 1]';
x2 = rho2*cos(theta2)*sin(pi/6);
y2 = rho2*sin(theta2)*sin(pi/6);
z2 = rho2*ones(size(theta2))*cos(pi/6) + 0.35;
cone = surf(x2,y2,z2);
axis square
c1 = [1 1]'*cos(theta2);
set(cone,'CData',c1, 'FaceColor', [0.2 0 0])
view(45,6)
hold off

end
