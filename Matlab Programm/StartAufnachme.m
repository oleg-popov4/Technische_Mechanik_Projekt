clear, close all

gitter = {'mesh==test==05'};


%Zeitangabe
ZeitEnde = 100;
CFL = 0.5;
Plotaktuallisirung = 0.1;

% resolution = [0 0 2560 1440];
resolution = [0 0 1920 1080];

aufnachme = 'j';

for itgitter = 1 : length(gitter)
    gitter_name = gitter{itgitter};
    run('Aufnachme.m');
end
