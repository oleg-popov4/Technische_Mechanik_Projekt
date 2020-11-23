function [ u_temp, Ausbreitungsgeschw_max ] = rhs(dt,u,length_t,...
    vi_vektor_AlleDreicke,DreieckZuKante,NormalenVektoren,...
    Tn_Matrix_Cell,Inverse_Tn_Matrix_Cell,HLL_temp,nullvektor,...
    Ausbreitungsgeschw_max,KantenLange,DreickFlache,Kante_Reihenfolge)


Ausbreitungsgeschw_tem = [0 0 0];
u_temp = nullvektor;

for i = 1:length_t
    sum = [0; 0; 0];
    for m = 1:3
        %Flussfunktion H---------------------------------------------------
        nachbar_dreick_temp = DreieckZuKante(i,m);
        %Ist es ein Dreieck am Rand?
        if (nachbar_dreick_temp == 0)
            HLL_temp{i,m} = Inverse_Tn_Matrix_Cell{i,m}*[0; 0.5*u(1,i)^2; 0];
            vi_n = vi_vektor_AlleDreicke(:,i)' *NormalenVektoren{i,m};
            Ausbreitungsgeschw_tem(1,m) = max(abs([vi_n-sqrt(u(1,i)) vi_n vi_n+sqrt(u(1,i))]));
        elseif (nachbar_dreick_temp < i)
            %HLL_temp{i,m} = - HLL_temp{nachbar_dreick_temp,Kante_Reihenfolge(i,m)};
            HLL_temp{i,m} = - HLL_temp{nachbar_dreick_temp,Kante_Reihenfolge(i,m)};
        else
            %Bestimung von u an der Kante von Dreieck i
            normal_vek_temp = NormalenVektoren{i,m};
            vi_n = vi_vektor_AlleDreicke(:,i)' *normal_vek_temp;
            vj_n = vi_vektor_AlleDreicke(:,nachbar_dreick_temp)' *normal_vek_temp;
            
            SL = min([vi_n-sqrt(u(1,i)); vj_n-sqrt(u(1,nachbar_dreick_temp))]);
            SR = max([vi_n+sqrt(u(1,i)); vj_n+sqrt(u(1,nachbar_dreick_temp))]);
            wi = Tn_Matrix_Cell{i,m}*u(:,i);
            wj = Tn_Matrix_Cell{i,m}*u(:,nachbar_dreick_temp);
            %[wj(2); wj(2)^2/wj(1)+0.5*wj(1)^2; wj(2)*wj(3)/wj(1)];
            %Hi = f( wi );
            
            Hi = [wi(2); wi(2)^2/wi(1)+0.5*wi(1)^2; wi(2)*wi(3)/wi(1)];
            %Hj = f( wj );
            Hj = [wj(2); wj(2)^2/wj(1)+0.5*wj(1)^2; wj(2)*wj(3)/wj(1)];
            if (SL >= 0)
                HLL_temp{i,m} = Inverse_Tn_Matrix_Cell{i,m}*Hi;
            elseif (SL <= 0 && SR >= 0)
                HLL_temp{i,m} = (SR*Hi - SL*Hj + SL*SR*(wj-wi))/(SR-SL);
                HLL_temp{i,m} = Inverse_Tn_Matrix_Cell{i,m}*HLL_temp{i,m};
            elseif (SR <= 0)
                HLL_temp{i,m} = Inverse_Tn_Matrix_Cell{i,m}*Hj;
            end
            Ausbreitungsgeschw_tem(1,m) = max(abs([vi_n-sqrt(u(1,i)) vi_n vi_n+sqrt(u(1,i))]));
        end
        %Ausbreitungsgeschw_max(i) = max(Ausbreitungsgeschw_tem);
        sum = sum + KantenLange(i,m)*HLL_temp{i,m};
    end
    Ausbreitungsgeschw_max(i) = max(Ausbreitungsgeschw_tem);
    %u_temp(:,i) = ( dt/DreickFlache(i) ) * KantenLange(i,:)*[HLL_temp{i,:}]';
    u_temp(:,i) = -( 1/DreickFlache(i) ) * sum;
end

end

