function [updicc] = patch_ups(dicc)
%PATCH_UPS sobremuestrear parches de un diccionario.
%Input:
%   dicc:diccionario en LR 27 X n_parches
%Output:
%   updicc:diccionario sobremuetreado 216 X n_parches

[~,n_parches]=size(dicc);
updicc=zeros(6^3,n_parches);

[Xq,Yq,Zq]=meshgrid((1:2/5:3),(1:2/5:3),(1:2/5:3));
for i=1:n_parches
    pch3d=reshape(dicc(:,i),[3,3,3]); %parche 3x3x3
    pchint=interp3(pch3d,Xq,Yq,Zq); %interpolando a parche 6x6x6
    updicc(:,i)=reshape(pchint,[6*6*6,1]);
    updicc(:,i)=updicc(:,i)/norm(updicc(:,i)); %vectores unitarios
end
end

