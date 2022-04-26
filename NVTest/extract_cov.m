function [vec] = extract_cov(Mx)
%EXTRACT_COV Extrae los valores de varianza y covarianza en un vector
%   Se descartan los valores repetidos de covarianza
%   Input:
%       Mx: Matriz cuadrada de covarianza (simétrica)
%   Output:
%       vec: Vector con valores de varianza y covarianza
[rows,cols]=size(Mx);
Mx_tri=triu(Mx);
vec_tri=reshape(Mx_tri,[rows*cols,1]);
ix=find(vec_tri);
vec=vec_tri(ix);
end

