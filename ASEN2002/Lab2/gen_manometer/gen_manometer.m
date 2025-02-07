%Josh Mellin
%This function takes the pressure in psi as an input and returns the x and
%y points needed for generating the manometer plots.  Also returns the
%difference in height of the fluid in meters

function [xp, yp, diff] = gen_manometer(pressure, rho)
%set up the initial radii and theta linspace
n = 50;
t = linspace(pi,2*pi,4*n);
g = 9.81; %force of gravity [m/s^2]
height = 0.25;
r1 = height/2;
r2 = height;

%get the base of the manometers
x1 = r1*cos(t);
y1 = r1*sin(t);
x2 = r2*cos(t);
y2 = r2*sin(t);

%generate the points for the walls
xlines = -r1*ones(1,n);
xlines = [xlines r1*ones(1,n)];
xlines = [xlines -r2*ones(1,n)];
xlines = [xlines r2*ones(1,n)];

ylines = linspace(0,height,n);
ylines = [ylines linspace(0,height,n)];
ylines = [ylines linspace(0,height,n)];
ylines = [ylines linspace(0,height,n)];

%create the points for the levels
diff = pressure / (rho*g);
level1 = height/2+diff;
level2 = height/2-diff;

%fill in the bottom circle first
np = 100;
r = linspace(r1,r2,np);
t1 = linspace(0,level1,np);
t2 = linspace(0,level2,np);
xin = [];
yin = [];
xt1 = []; %left side
yt1 = []; %left side
xt2 = []; %right side
yt2 = []; %right side
for k = 1:length(r)
    xin = [xin r(k)*cos(t)];
    yin = [yin r(k)*sin(t)];
    xt1 = [xt1 linspace(-r2, -r1, n)];
    xt2 = [xt2 linspace(r1, r2, n)];
    yt1 = [yt1 t1(k)*ones(1,n)];
    yt2 = [yt2 t2(k)*ones(1,n)];
end


%concatenate all points into one vector
xp = [x1 x2 xlines xin xt1 xt2];
yp = [y1 y2 ylines yin yt1 yt2];

yp = yp*1000; %convert from m to mm


end


