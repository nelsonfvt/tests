function [B_sim] = pch_sim(dd,DD)
%PCH_SIM Función para encontrar parches semejantes
%Input:
%   dd: Conjunto de vectores patrón
%   DD: Conjunto total de vectores
%Output:
%   B_sim: Matriz
%

[drows,dcols]=size(dd);
[~,Dcols]=size(DD);

dists=zeros(dcols,Dcols);
B_sim=zeros(drows,dcols,dcols);
for i=1:dcols
    for j=1:Dcols
        dists(i,j)=dot(dd(:,i),DD(:,j));
    end
    [~,I]=sort(dists(i,:),'descend');% ordenando desde el más parecido
    for k=1:dcols
        B_sim(:,i,k)=DD(:,I(k));% llenando arreglo con parches más parecidos
    end
end

end

