%% Shallow_water_test mesh generation
xc_1 = 12;
yc_1 = 2.5;
r_1 = 0.5*pi;
fd_Tepm = @(p) dcircle(p,xc_1,yc_1,r_1);

xc_2 = 3.5;
yc_2 = 0;
r_2 = 0.5*pi;
fd_Tepm_2 = @(p) dcircle(p,xc_2,yc_2,r_2);

xc_3 = 7.5;
yc_3 = 5;
r_3 = 0.5*pi;
fd_Tepm_3 = @(p) dcircle(p,xc_3,yc_3,r_3);

%fd = @(p)
%dunion(drectangle(p,-5,0,0,5),ddiff(drectangle(p,0,15,0,5),fd_Tepm(p)));dintersect
%ddiff(drectangle(p,0,15,0,5),dintersect(drectangle(p,0,15,0,5),fd_Tepm_2(p)))
fd = @(p) dunion(drectangle(p,-5,0,0,5), ddiff( dintersect(ddiff(drectangle(p,0,15,0,5),dintersect(drectangle(p,0,15,0,5),fd_Tepm_3(p))),ddiff(drectangle(p,0,15,0,5),dintersect(drectangle(p,0,15,0,5),fd_Tepm_2(p))) ),fd_Tepm(p)) );
pfix = [0,0; 15,0; 15,5; 0,5; 2, 1.5; 4,1.5; 4,3.5; 2,3.5];
pfix = [pfix; -5,0; -5,5];

% Decrease dx to refine mesh
dx = .1;

pfix = [pfix; [zeros(size((0:dx:5)')) (0:dx:5)']];
[p,t]=distmesh2d(fd,@huniform,dx,[-5,0;15,5],unique(pfix,'rows'));
%set(gca,'visible','on')
% Uncomment to save mesh in file mesh.mat
save spektakulearer t p