function [vals] = covars(dd,DD)
%COVARS Función para calcular las covarianzas entre parches base y otros
%Input:
%   cov_base: parches base
%   DD: conjunto de vectores para comparar
%Output:
%   vals: valores de las covarianzas
cov_base = extract_cov(cov(dd));
%cov_base=cov_base/norm(cov_base);
[~,~,n_per]=size(DD);
vals=zeros(1,n_per);

for k=1:n_per
    covs=extract_cov(cov(DD(:,:,k)));
    %covs=covs/norm(covs);
    vals(k)=norm(cov_base-covs);
end

end

