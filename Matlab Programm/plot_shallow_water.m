function [VideoFrame] = plot_shallow_water(u,time,t,p,phimin, phimax, vmin, vmax,gitter_name,hfig)
% PLOT_SHALLOW_WATER Plot geopotential and norm of velocity.
%      u: Matrix of size 3 x length(t) which contains the conservative
%         variables of the shallow water equations for all triangles
%   time: Current time
%      t: Triangles of the mesh
%      p: Coordinates of the mesh vertices
% phimin: Lower limit of the z axes when phi is plotted
% phimax: Upper limit of the z axes when phi is plotted
%   vmin: Lower limit of the z axes when v is plotted
%   vmax: Upper limit of the z axes when v is plotted

subplot(221)
trisurf_pc(t,p(:,1),p(:,2),u(1,:)');
set(gca,'CLim',[phimin phimax]);
shading flat
axis equal
axis tight
zlim([phimin phimax])
title([sprintf('Time: t = %.2f\n',time) gitter_name])
ylabel('\phi','Rotation',0,'FontSize',22,'FontWeight','bold');
%colorbar

subplot(223)
trisurf_pc(t,p(:,1),p(:,2),sqrt((u(2,:)'./u(1,:)').^2+(u(3,:)'./u(1,:)').^2));
shading flat
axis equal
axis tight
set(gca,'CLim',[vmin, vmax]);
zlim([vmin vmax])
title([sprintf('Time: t = %.2f\n',time) gitter_name])
ylabel('||v||','Rotation',0,'FontSize',22,'FontWeight','bold'); 
colormap(jet)
%colorbar

subplot(222)
trisurf_pc(t,p(:,1),p(:,2),u(1,:)');
set(gca,'CLim',[phimin phimax]);
shading flat
axis equal
axis tight
zlim([phimin phimax])
view(2)
title([sprintf('Time: t = %.2f\n',time) gitter_name])
y= ylabel('\phi','Rotation',0,'FontSize',22,'FontWeight','bold');
set(y, 'position', get(y,'position')-[1,0,0]); 
%colorbar

subplot(224)
trisurf_pc(t,p(:,1),p(:,2),sqrt((u(2,:)'./u(1,:)').^2+(u(3,:)'./u(1,:)').^2));
shading flat
axis equal
axis tight
set(gca,'CLim',[vmin, vmax]);
zlim([vmin vmax])
view(2)
title([sprintf('Time: t = %.2f\n',time) gitter_name])
y= ylabel('||v||','Rotation',0,'FontSize',22,'FontWeight','bold');
set(y, 'position', get(y,'position')-[1,0,0]); 
colormap(jet)
%colorbar




try
    im = print(hfig,'-RGBImage','-r0');
    VideoFrame = im2frame(im);
catch Me
    VideoFrame = 0;
end


