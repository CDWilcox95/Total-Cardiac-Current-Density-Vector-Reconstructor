function setupfigure

load colormap.mat
% CMap='jet';
shading interp ; % make it more smooth
colormap(CMap);

axis xy;
axis normal;
axis off;
end

