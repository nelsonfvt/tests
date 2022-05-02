function [LRdicc] = patch_dws(HRdicc)
%PATCH_DWS baja la resolución de un diccionario de 6*6*6 a 3*3*3
%Input
%   HRdicc:diccionario en HR 216 X n_patch
%Output
%   LRDicc:diccionario en LR 27 X n_patch

[~,n_patch]=size(HRdicc);
LRdicc=zeros(3^3,n_patch);

[Xq,Yq,Zq]=meshgrid([1 3.5 6],[1 3.5 6],[1 3.5 6]);
for i=1:n_patch
    pch3d=reshape(HRdicc(:,i),[6,6,6]); %parche 6x6x6
    pchint=interp3(pch3d,Xq,Yq,Zq); %interpolando a parche 3x3x3
    LRdicc(:,i)=reshape(pchint,[3*3*3,1]);
    LRdicc(:,i)=LRdicc(:,i)/norm(LRdicc(:,i)); %vectores unitarios
end
end

