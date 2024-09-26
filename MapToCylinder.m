function h=MapToCylinder(R, Vq,bg)


load("SurfaceMap_Vertices.mat");    load("SurfaceMap_Faces.mat");
% 
% % Visualize cylinder
% figure('color','w')
% axis equal
% h=patch('Faces',F,'Vertices',V);
% set(h,'FaceColor',0.75*[1 1 1],'FaceAlpha',0.8,'EdgeColor','k')
% view([20 20])
% xlabel('X','FontSize',20,'FontWeight','bold')
% ylabel('Y','FontSize',20,'FontWeight','bold')
% zlabel('Z','FontSize',20,'FontWeight','bold')

L=16;
theta_l=(2*pi/L).*[1:L]'; Z(1:L)=0.04;
elec_pts=[R.*cos(theta_l), R.*sin(theta_l),0.04.*ones(L,1); R.*cos(theta_l), R.*sin(theta_l), 0.11.*ones(L,1)];

% Visualize G
% axis equal
Vq(isnan(Vq))=bg;

V(:,3)=(10^-2).*V(:,3);
h=patch('Faces',F,'Vertices',V,'FaceVertexCData',Vq(:),'FaceColor','interp');
set(h,'FaceAlpha',0.95,'EdgeColor','none');
hold on;
view([30 40])

shift=0.005;
% plot3(R.*cos(elec_pts(:,1)), R.*sin(elec_pts(:,1)), elec_pts(:,2)+shift, 'ko', 'MarkerSize', 15, 'MarkerFaceColor','w');
plot3(elec_pts(:,1), elec_pts(:,2), elec_pts(:,3), 'ko', 'MarkerSize', 15, 'MarkerFaceColor','w');

for i=1:length(elec_pts)
    text(elec_pts(i,1)-shift, elec_pts(i,2)-shift, elec_pts(i,3), num2str(i));
end
hold on;
% xlabel('X','FontSize',20,'FontWeight','bold')
% ylabel('Y','FontSize',20,'FontWeight','bold')
% zlabel('Z','FontSize',20,'FontWeight','bold')
% camlight('headlight'), lighting phong
% frame=getframe(f);
% writeVideo(vidfile, frame);
end