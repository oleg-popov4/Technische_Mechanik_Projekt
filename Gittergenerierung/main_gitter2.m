% Shallow_water_test mesh generation
clc
clear all
close all

%Beim 'ja' - wird gespeichert sonst nicht
speichern = 'j';
Grid_name = 'Rechteck==im==rechteck';

% Decrease dx to refine mesh
dx_vektor = [0.025 0.01];


for iter = 1:length(dx_vektor)
    grid_name = Grid_name;
    dx = dx_vektor(iter);
    
    %Fuer Beschlenigung  beim distmesh2d zeile 92 komentieren
    %x1,y1,x2,y2
    r1 = [-5,1,0,6];
    
    %r2 ist ergänzung an r1!!!
    r2 = [0,0,10,7];
    
    %rausgeschnitene Rechteck
    r3 = [5,1,6,3];
    
    %=drectangle(p,x1,x2,y1,y2)
    %dunion
    %ddiff
    fd = @(p) ...
    ddiff ( ...
            dunion(drectangle(p,r1(1), r1(3),r1(2),r1(4)),drectangle(p,r2(1), r2(3),r2(2),r2(4) ) ),...
            drectangle(p,r3(1), r3(3),r3(2),r3(4)) ...
           );
    
    %Fixire Punkte am Rand fuer Rechteck 2
    pfix =   [r2(1),r2(2); r2(3),r2(2); r2(3),r2(4); r2(1),r2(4)];
    %Fixire Punkte am Rand fuer Rechteck 3
    pfix = [pfix; r3(1),r3(2); r3(3),r3(2); r3(3),r3(4); r3(1),r3(4)];
    
    
    %Fixire Punkte am Rand fuer Rechteck 1
    pfix = [pfix; r1(1), r1(2); r1(1),r1(end)];
    
    
    %Setze zusätzliche Punkte auf verbindungsstrecke von r1 und r2
    [x_strecke,fx_strecke]= gerade_funkt([r1(3),r1(2)],[r1(3),r1(end)],dx);
    pfix = [pfix; [x_strecke' fx_strecke']];
    
    
    %[p,t]=distmesh2d(fd,@huniform,dx,[r1(1),r1(2);r2(3),r2(4)],unique(pfix,'rows'));
    [p,t]=distmesh2d(fd,@huniform,dx,[-100,-100;100,100],unique(pfix,'rows'));
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
