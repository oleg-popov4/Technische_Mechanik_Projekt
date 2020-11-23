% Shallow_water_test mesh generation
clc
clear all
close all

%Beim 'ja' - wird gespeichert sonst nicht
speichern = 'j';
Grid_name = 'mesh_test_01';

% Decrease dx to refine mesh
dx_vektor = [.5];


for iter = 1:length(dx_vektor)
    grid_name = Grid_name;
    dx = dx_vektor(iter);
    
    %Fuer Beschlenigung  beim distmesh2d zeile 92 komentieren
    %x1,y1,x2,y2
    r1 = [-1,0,0,6];
    
    %r2 ist ergänzung an r1!!!
    r2 = [0,0,6,6];
    
    %rausgeschnitene Rechteck
    r3 = [3,1,5,1.5];
    
    %Hinzugefugte Kreis
    kreis1_x = 6;
    kreis1_y = 3;
    kreis1_radius = 3;
    
    %Hinzugefugte Kreis 2
    kreis2_x = -1;
    kreis2_y = 3;
    kreis2_radius = 3;
    
    %Rausgeschnitene Kreis1
    kreis3_x = 2.5;
    kreis3_y = 3;
    kreis3_radius = 0.5;
    
    %Rausgeschnitene Kreis2
    kreis4_x = 5.5;
    kreis4_y = 3;
    kreis4_radius = 0.5;
    
    %=drectangle(p,x1,x2,y1,y2)
    %dunion
    %ddiff
    fd = @(p) ...
    ddiff ( ...
    ddiff ( ...
    ddiff ( ...
    dunion( ...
            dunion(drectangle(p,r1(1), r1(3),r1(2),r1(4)),...
                    dunion(...
                           drectangle(p,r2(1), r2(3),r2(2),r2(4)),...
                           dcircle(p, kreis1_x, kreis1_y,kreis1_radius)   ) ),...
            dcircle(p, kreis2_x, kreis2_y,kreis2_radius) ),...
            drectangle(p,r3(1), r3(3),r3(2),r3(4))       ), ...
    dcircle(p, kreis3_x, kreis3_y,kreis3_radius) ), ...
    dcircle(p, kreis4_x, kreis4_y,kreis4_radius) );
    
    %Fixire Punkte am Rand fuer Rechteck 2
    pfix =   [r2(1),r2(2); r2(3),r2(2); r2(3),r2(4); r2(1),r2(4)];
    %Fixire Punkte am Rand fuer Dreick 3
    pfix = [pfix; r3(1),r3(2); r3(3),r3(2); r3(3),r3(4); r3(1),r3(4)];
    
    
    %Fixire Punkte am Rand fuer Rechteck 1
    pfix = [pfix; r1(1), r1(2); r1(1),r1(end)];
    
    
    %Setze zusätzliche Punkte auf verbindungsstrecke von r1 und r2
    [x_strecke,fx_strecke]= gerade_funkt([r1(2),r1(3)],[r1(2),r1(end)],dx);
    pfix = [pfix; [x_strecke' fx_strecke']];
    
    
    %[p,t]=distmesh2d(fd,@huniform,dx,[r1(1),r1(2);r2(3),r2(4)],unique(pfix,'rows'));
    [p,t]=distmesh2d(fd,@huniform,dx,[-10,-10;10,10],unique(pfix,'rows'));
    set(gca,'visible','on')
    
    
    
    
    
    
    grid_name = [grid_name, '_{' , num2str(dx),'}==nodes=',num2str(size(p,1)),...
        '==elements=',num2str(size(t,1)),...
        '==quality=',num2str(min(simpqual(p,t)))];
    grid_name = strrep(grid_name,'.',',');
    grid_name = [grid_name '.mat'];
    % Uncomment to save mesh in file mesh.mat
    if ( strcmp(speichern,'ja') )
        save(grid_name,'t','p');... grid_name t p
    end

end


function [x_vektor, fvektor ] = gerade_funkt ( p0, p1, dx )
if ( p0(1) ~= p1(1) )
    I(1) = min(p0(1),p1(1));
    I(2) = max(p0(1),p1(1));
    x_vektor = I(1):dx:I(2);
    nenner = p0(1) - p1(1);
    a = -(-p0(2) + p1(2))/nenner;
    b = -(p0(2)*p1(1) - p0(1)*p1(2) )/nenner;
    b_vektor = b*ones(size(x_vektor));
    fvektor = a*x_vektor + b_vektor;
else
    I(1) = min(p0(2),p1(2));
    I(2) = max(p0(2),p1(2));
    fvektor = I(1):dx:I(2);
    x_vektor = p0(1)*ones(size(fvektor));
end

end







% fd = @(p) dunion(drectangle(p,r1(1), r1(3),r1(2),r1(4)),...
%      dunion(...
%     ddiff(drectangle(p,r2(1), r2(3),r2(2),r2(4)),...
%     drectangle(p,r3(1), r3(3),r3(2),r3(4)) ), ... , ...
%      dcircle(p,6,3,3)         )...
%     );
%
%Komische schwimmbecken
% fd = @(p) ...
%   dunion( ...
%     dunion(drectangle(p,r1(1), r1(3),r1(2),r1(4)),...
%        dunion(...
%               ddiff(drectangle(p,r2(1), r2(3),r2(2),r2(4)),...
%                      drectangle(p,r3(1), r3(3),r3(2),r3(4)) ), ... , ...
%               dcircle(p, kreis1_x, kreis1_y,kreis1_radius)         ) ),...
%      dcircle(p, kreis2_x, kreis2_y,kreis2_radius)...
%           );




% [x_var, fx_var] = gerade_funkt([r3(1),r3(2)],[r3(1),r3(end)],dx);
% pfix = [pfix; [x_var' fx_var']];
%
% [x_var, fx_var] = gerade_funkt([r3(1),r3(end)],[r3(3),r3(4)],dx);
% x_var = x_var(2:end-1);
% fx_var = fx_var(2:end-1);
% pfix = [pfix; [x_var' fx_var']];
%
% [x_var, fx_var] = gerade_funkt([r3(3),r3(4)],[r3(3),r3(2)],dx);
% pfix = [pfix; [x_var' fx_var']];
%
% [x_var, fx_var] = gerade_funkt([r3(1),r3(2)],[r3(3),r3(2)],dx);
% x_var = x_var(2:end-1);
% fx_var = fx_var(2:end-1);
% pfix = [pfix; [x_var' fx_var']];

