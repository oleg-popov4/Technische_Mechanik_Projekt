clear, close all


%Zeitangabe
ZeitEnde = 1;
CFL = 0.5;


%Gitter laden
gitter_name = 'gitter=smailli_{0,1}==nodes=7858==elements=15276==quality=0,73495';
load([gitter_name '.mat'])
gitter_name = strrep(gitter_name,'==',' ');

%aufnachme - ja speichert sonst nicht
aufnachme = 'j';
Plotaktuallisirung = 0.1;

resolution = [0 0 1920 1080];
% Create a new figure
hfig = figure('units','pixels','position',resolution,'Visible','on', ... % New figure is invisible
    'PaperPositionMode','auto');

if ( strcmp(aufnachme,'ja')  )
    % Create a new VideoWriter object (an empty video file). Use whatever format you want,
    % but in my experience MP4 (with H.264 codec) is by far the best. Please stop using AVI.
    hvid = VideoWriter([gitter_name '.avi']);
    %     vidObj = VideoWriter([gitter_name '.avi']);
    %     vidObj.FrameRate = 10;
    %     writerObj.Height=1080;
    %     writerObj.Width=1920;
    %     open(vidObj);
    % Full quality, because why not?
    set(hvid,'Quality',100);
    
    % Set the frame rate
    set(hvid,'FrameRate',10);
    % Open the object for writing
    open(hvid);
    framepar.resolution = [resolution(3) resolution(4)];
end

%Ausgabe der Information uber Elemente
fprintf('%d nodes, %d elements, min quality %.2f\n',size(p,1),size(t,1),min(simpqual(p,t)));

zeitmessung = tic;


AnzahlPunkte = size(t);
DreieckKoordinaten = zeros(2,AnzahlPunkte(2));
u0 = zeros(3,AnzahlPunkte(1));
DreickFlache = u0;
KantenLange = zeros(AnzahlPunkte);
NormalenVektoren =  cell(AnzahlPunkte);
DreieckZuKante = KantenLange;
%Kante_Reihenfolge = sparse(AnzahlPunkte(1),AnzahlPunkte(2));
Kante_Reihenfolge = zeros(AnzahlPunkte(1),AnzahlPunkte(2));
k_Abstand_Schwerpunkts = zeros(AnzahlPunkte(1),1);

%-----------------------------------------------------

Tn_Matrix_Cell =  cell(AnzahlPunkte);
Inverse_Tn_Matrix_Cell =  cell(AnzahlPunkte);
Tn_Matrix = @(n) [1 0 0; 0 n(1) n(2); 0 -n(2) n(1)];
Tn_Matrix_Inverse1 = @(n) [1 0 0; 0 n(1) -n(2); 0 n(2) n(1)];
f = @(u) [u(2); u(2)^2/u(1)+0.5*u(1)^2; u(2)*u(3)/u(1)];
Ausbreitungsgeschw_max = zeros(AnzahlPunkte(1),1);

%-------------------------------------------------------

%Alle Dreiecke durchgehen
for i= 1:length(t)
    %Nimmt Knoten des Dreiecks für Bestimmung der Koordinaten der Punkte
    %P1, P2, P3 vom Dreieck.
    DreieckKnoten = t(i,:);
    %Berechnnung der Dreieck Koordinaten
    for j = 1:length(DreieckKnoten)
        DreieckKoordinaten(:,j) = p(DreieckKnoten(j),:);
    end
    %Definition der Anfangsdaten
    [u0(:,i), DreickFlache(i), zustand_links, zustand_rechts] = Mittelwerte(DreieckKoordinaten);
    k_Abstand_Schwerpunkts(i) = comp_k(DreieckKoordinaten');
    %Kanten werden in fest von mir vorgeschriebene Reihenfolge bearbeitet
    %Zuerst P1-P2, P1-P3, P2-P3
    P1 = DreieckKoordinaten(:,1);
    P2 = DreieckKoordinaten(:,2);
    P3 = DreieckKoordinaten(:,3);
    %--------------
    KantenLange(i,:)= [norm(P1-P2,2) norm(P1-P3,2) norm(P2-P3,2)];
    %--------------
    %Orthogonalprojektion Quelle
    %https://de.wikipedia.org/wiki/Orthogonalprojektion
    
    %P1-P2
    u = P2-P1;
    Pg = P1 + ( (P3-P1)'*u/(u'*u) )*u;
    NormalenVektoren{i,1} = (Pg-P3)/norm(Pg-P3);
    
    Tn_Matrix_Cell{i,1} = Tn_Matrix(NormalenVektoren{i,1});
    Inverse_Tn_Matrix_Cell{i,1} = Tn_Matrix_Inverse1(NormalenVektoren{i,1});
    %     if(abs(u'*NormalenVektoren{i,1})>=1e-15)
    %         warning(['P1-P2, bei dreieck = ', num2str(i)]);
    %         u'*NormalenVektoren{i,1}
    %         pause;
    %     end
    
    %P1-P3
    u = P3-P1;
    Pg = P1 + ( (P2-P1)'*u/(u'*u) )*u;
    NormalenVektoren{i,2} = (Pg-P2)/norm(Pg-P2);
    
    Tn_Matrix_Cell{i,2} = Tn_Matrix(NormalenVektoren{i,2});
    Inverse_Tn_Matrix_Cell{i,2} = Tn_Matrix_Inverse1(NormalenVektoren{i,2});
    %     if(abs(u'*NormalenVektoren{i,2})>=1e-15)
    %         warning(['P1-P3, bei dreieck = ', num2str(i)]);
    %         pause;
    %     end
    
    %P2-P3
    u = P3-P2;
    Pg = P2 + ( (P1-P2)'*u/(u'*u) )*u;
    NormalenVektoren{i,3} = (Pg-P1)/norm(Pg-P1);
    
    Tn_Matrix_Cell{i,3} = Tn_Matrix(NormalenVektoren{i,3});
    Inverse_Tn_Matrix_Cell{i,3} = Tn_Matrix_Inverse1(NormalenVektoren{i,3});
    %     if(abs(u'*NormalenVektoren{i,3})>=1e-15)
    %         warning(['P2-P3, bei dreieck = ', num2str(i)]);
    %         pause;
    %     end
    
    v_temp = [u0(2,i)./u0(1,i) u0(3,i)./u0(1,i)];
    
    lambda2_1 = v_temp*NormalenVektoren{i,1};
    lambda1_1 = lambda2_1-sqrt(u0(1,i));
    lambda3_1 = lambda2_1+sqrt(u0(1,i));
    
    lambda2_2 = v_temp*NormalenVektoren{i,2};
    lambda1_2 = lambda2_1-sqrt(u0(1,i));
    lambda3_2 = lambda2_1+sqrt(u0(1,i));
    
    lambda2_3 = v_temp*NormalenVektoren{i,3};
    lambda1_3 = lambda2_1-sqrt(u0(1,i));
    lambda3_3 = lambda2_1+sqrt(u0(1,i));
    
    Ausbreitungsgeschw_max(i) = max(abs( [lambda2_1 lambda1_1 lambda3_1...
        lambda2_2 lambda1_2 lambda3_2 lambda2_3 lambda1_3 lambda3_3]) );
end

%Definition vom Triangulirung
TR = triangulation(t,p(:,1),p(:,2),0*p(:,1));

%Plot der Anfangsbedingung

%Plot der Randwerte--------------------------------------------------------
phimin = min(zustand_links, zustand_rechts)-1;
phimax = max(zustand_links, zustand_rechts)+1;
vmin = -1;
vmax = max(zustand_links, zustand_rechts)*0.5;

if ( strcmp(aufnachme,'ja')  )
    figure(hfig)
end
F = plot_shallow_water(u0,0,t,p,phimin,phimax,vmin,vmax,gitter_name,hfig);
if ( strcmp(aufnachme,'ja')  )
    % Convert the figure to a video frame.
    % The built-in function for this is GETFRAME, which has a variety of annoying features.
    % fig2frame is a drop-in replacement for getframe that avoids most (all?) the annoyance.
    
    %F = fig2frame(hfig,framepar); % <-- Use this
    % F = getframe(hfig); % <-- Not this.
    
    % Add the frame to the video object
    writeVideo(hvid,F);
    savefig([gitter_name 'Zeit=0.fig'])
end

%Finite-Volumen-Verfahren---------------------------------------------------

%Gehe alle Dreieck durch, erzeuge eine Lise die durgelaufen wird
DreieckeAnKante = cell(3,1);
for i = 1:length(t)
    
    DreieckKnoten = t(i,:);
    DreieckeAnKante{1} = edgeAttachments(TR,DreieckKnoten(1),DreieckKnoten(2));
    DreieckeAnKante{2} = edgeAttachments(TR,DreieckKnoten(1),DreieckKnoten(3));
    DreieckeAnKante{3} = edgeAttachments(TR,DreieckKnoten(2),DreieckKnoten(3));
    for m = 1:3
        if length( DreieckeAnKante{m}{1} )==1
            %Kante liegt am Rand der Triangulirung
            DreieckZuKante(i,m) = 0;
        else
            %Ein Dreick ist i, ich brauche aber anderen Dreieck.
            var1 = DreieckeAnKante{m}{1};
            if( var1(1) ~= i )
                DreieckZuKante(i,m) = var1(1);
                if var1(1) < var1(2)
                    Kante_Reihenfolge(i,m) = find(DreieckZuKante(var1(1),:)==var1(2));
                end
            else
                DreieckZuKante(i,m) = var1(2);
                if var1(2) < var1(1)
                    Kante_Reihenfolge(i,m) = find(DreieckZuKante(var1(2),:)==var1(1));
                end
            end
        end
    end
end

%Eigentliche Verfahren Summe über alle Dreiecke.-------------------
u = u0;

%Zeitschleife----------------------------------------------------
dt = min(k_Abstand_Schwerpunkts./Ausbreitungsgeschw_max);
if dt > Plotaktuallisirung
    dt = Plotaktuallisirung;
    Plot_laufvariable = 2;
else
    Plot_laufvariable = 1;
end
nullvektor = 0*u0;
nachbar_dreick_temp = 0;
HLL_temp = cell(AnzahlPunkte);
Zeit_temp = dt;
iter = 2;
Alle_Zeitschritte = zeros(100,1);
Alle_Zeitschritte(1) = dt;
while(Zeit_temp <= ZeitEnde)
    
    vi_vektor_AlleDreicke = [u(2,:)./u(1,:); u(3,:)./u(1,:)];
    length_t = length(t);
    
    [u_temp, Ausbreitungsgeschw_max] = rhs_mex(dt, u,length_t,...
        vi_vektor_AlleDreicke,DreieckZuKante,NormalenVektoren,...
        Tn_Matrix_Cell,Inverse_Tn_Matrix_Cell,HLL_temp,nullvektor,...
        Ausbreitungsgeschw_max,KantenLange,DreickFlache,Kante_Reihenfolge );
    
    %Berechnung der gesuchten u
    %Explizite Euler
    u = u + dt * u_temp;
    
    
    
    %Abfrage ob die u geplottet werden soll
    if Zeit_temp >= Plotaktuallisirung*Plot_laufvariable
        drawnow
        F = plot_shallow_water(u,Zeit_temp,t,p,phimin,phimax,vmin,vmax,gitter_name,hfig);
        if ( strcmp(aufnachme,'ja')  )
            % Convert the figure to a video frame.
            % The built-in function for this is GETFRAME, which has a variety of annoying features.
            % fig2frame is a drop-in replacement for getframe that avoids most (all?) the annoyance.
            
            %F = fig2frame(hfig,framepar); % <-- Use this
            % F = getframe(hfig); % <-- Not this.
            
            % Add the frame to the video object
            writeVideo(hvid,F);
            if(max(Zeit_temp==0:ZeitEnde))
                savefig([gitter_name 'Zeit=' num2str(Zeit_temp) '.fig'])
            end
        end
        Plot_laufvariable = Plot_laufvariable+1;
    end
    
    %Berechnung der neue Zeitschrittweitte
    dt_temp = min(k_Abstand_Schwerpunkts./Ausbreitungsgeschw_max);
    %Aktuallisieren der berechnete Zeit
    Zeit_temp = Zeit_temp + dt_temp;
    %Anpasung der Zeitschrittweitte an Plott ausgabe
    if(Zeit_temp >= Plotaktuallisirung*Plot_laufvariable)
        Zeit_temp = Zeit_temp - dt_temp;
        dt_temp = Plotaktuallisirung*Plot_laufvariable - Zeit_temp;
        Zeit_temp = Zeit_temp + dt_temp;
    end
    
    %Anpassung der Zeit an Zeit_Ende. Damit nicht uberschrietten wird
    if (Zeit_temp > ZeitEnde)
        if ((Zeit_temp - dt_temp)~=ZeitEnde)
            Zeit_temp = Zeit_temp - dt_temp;
            dt_temp = ZeitEnde - Zeit_temp;
            if dt_temp > 0
                Zeit_temp = Zeit_temp + dt_temp;
            end
        end
    end
    
    Alle_Zeitschritte(iter)=dt;
    iter = iter+1;
    dt = dt_temp;
end


if ( strcmp(aufnachme,'ja')  )
savefig([gitter_name 'Zeit=' num2str(Zeit_temp) '.fig'])
% Close the video object. This is important! The file may not play properly if you don't close it.
close(hvid);
close all
end


tElapsed = toc(zeitmessung)

% figure
% Alle_Zeitschritte = Alle_Zeitschritte(1:iter-1);
% plot(1:length(Alle_Zeitschritte),Alle_Zeitschritte,'.')
% title('Zeitschrittweite')
% xlabel('Iteration')
% ylabel('dt')










