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
n_sim = dcols/2;

dists=zeros(dcols,Dcols);
B_sim=zeros(drows,dcols,n_sim);
for i=1:dcols
    for j=1:Dcols
        %dists(i,j)=dot(dd(:,i),DD(:,j));
        dists(i,j)=norm(dd(:,i)-DD(:,j));
    end
    %[~,I]=sort(dists(i,:),'descend');% ordenando desde el más parecido
    [~,I]=sort(dists(i,:),'ascend');% ordenando desde el más parecido
    for k=1:n_sim
        B_sim(:,i,k)=DD(:,I(k));% llenando arreglo con atomos más parecidos
    end
end

end

