function [dicc_perm] = permutaciones(dicc)
%PERMUTACIONES Genera las permutaciones de los elementos del diccionario
%Input
%   dicc: diccionario o cleccion de parches (arreglo 3D)
%Output
%   dicc_perm: dicionario permutado

[rows,cols,n_elem]=size(dicc);

v=[1:n_elem]; % numero de elementos para permutar
C=perms(v);
C=[C sortrows(C,[1 2])];
% ns=1:20:length(C(:,1));
% C=C(ns,:);
n_per=length(C(:,1));
B_perm=zeros(rows,cols,n_per);

for i=1:n_per
    for j=1:cols
        vec=dicc(:,j,C(i,j));
        B_perm(:,j,i)=vec;
    end
end
dicc_perm=cat(3,B_perm,dicc);
end

