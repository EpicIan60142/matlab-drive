
load('FreeFall.mat', 'outputStateMatrix','timeVector');
x = outputStateMatrix(:,3)';
y = outputStateMatrix(:,4)';
z = zeros(size(x));

col = outputStateMatrix(:,2)'; %v-velocity as color
surface([x;x],[y;y],[z;z],[col;col],'facecol','no','edgecol','interp','linew',2);
xlabel('X position'); ylabel('Y position'); title('Y-velocity as color'); 
colorbar; colormap(jet); ax=gca; ax.FontSize= 16; ax.LineWidth = 1;