function [ ausgabe ] = u0( x, zustand_links, zustand_rechts )

if(x(1)<0)
    ausgabe = [zustand_links; 0; 0];
elseif(x(1)>0)
    ausgabe = [zustand_rechts; 0; 0];
else
    ausgabe = [0; 0; 0];
end

end

