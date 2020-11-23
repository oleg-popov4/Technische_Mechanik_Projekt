function [ Mittelwert, Flache, zustand_links, zustand_rechts] = Mittelwerte( Punkte )

%Zustand 0 und negative ist nicht erlaubt
zustand_links = 10;
zustand_rechts = 0.1;

%Punkt 1
x1 = Punkte (1,1);
y1 = Punkte (2,1);
%Punkt 2
x2 = Punkte (1,2);
y2 = Punkte (2,2);
%Punkt 3
x3 = Punkte (1,3);
y3 = Punkte (2,3);
%Abbildung Tau
Tau = @(xi,n) [x2-x1 x3-x1; y2-y1 y3-y1]*[xi;n] +[x1;y1];
%Jacobi Matrix
Jacobi = [x2-x1 x3-x1; y2-y1 y3-y1];
%Quadraturformel

Flache = abs(det(Jacobi))/2;

Integration =abs(det(Jacobi))* ( u0( Tau(0,0.5),zustand_links, zustand_rechts ) + u0( Tau(0.5,0),zustand_links, zustand_rechts )...
                + u0( Tau(0.5,0.5),zustand_links, zustand_rechts ) )/6;




Mittelwert = Integration/abs(Flache);

%Anpasung der Unstetigkeitsstelle
if (Mittelwert(1) ~= max(zustand_links,zustand_rechts) && Mittelwert(1) ~= min(zustand_links,zustand_rechts))
    if Mittelwert(1) > ( min(zustand_links,zustand_rechts) + 0.1)
        Mittelwert(1) = zustand_links;
    else
        Mittelwert(1) = zustand_rechts;
    end
end






end

