function display_cuboid(x,y,z)



dist_x = abs(x(end)-x(1));
dist_y = abs(y(end)-y(1));
dist_z = abs(z(end)-z(1));

[X1, Y1, Z1] = meshgrid(x([1 end]),y,z);
X1 = permute(X1,[2 1 3]); Y1 = permute(Y1,[2 1 3]); Z1 = permute(Z1,[2 1 3]);
X1(end+1,:,:) = NaN; Y1(end+1,:,:) = NaN; Z1(end+1,:,:) = NaN;
[X2, Y2, Z2] = meshgrid(x,y([1 end]),z);
X2(end+1,:,:) = NaN; Y2(end+1,:,:) = NaN; Z2(end+1,:,:) = NaN;
[X3, Y3, Z3] = meshgrid(x,y,z([1 end]));
X3 = permute(X3,[3 1 2]); Y3 = permute(Y3,[3 1 2]); Z3 = permute(Z3,[3 1 2]);
X3(end+1,:,:) = NaN; Y3(end+1,:,:) = NaN; Z3(end+1,:,:) = NaN;

figure('Renderer','opengl')
h = line([X1(:);X2(:);X3(:)], [Y1(:);Y2(:);Y3(:)], [Z1(:);Z2(:);Z3(:)]);
set(h, 'Color',[0.5 0.5 1], 'LineWidth',1, 'LineStyle','-.')

set(gca, 'Box','off', 'LineWidth',2, 'XTick',x, 'YTick',y, 'ZTick',z, ...
  'XLim',[x(1)-0.1*dist_x x(end)+0.1*dist_x], ...
  'YLim',[y(1)-0.1*dist_y y(end)+0.1*dist_y], ...
  'ZLim',[z(1)-0.1*dist_z z(end)+0.1*dist_z])
xlabel x, ylabel y, zlabel z
axis on
grid on
view(3), axis vis3d
camproj perspective, rotate3d on

end