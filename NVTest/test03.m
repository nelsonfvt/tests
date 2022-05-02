%% Inicio
[rows,cols]=size(DDHr);
t_vecs=find(WNueLi);

%% Upsampling parches LR
ddLR=DDLR(:,t_vecs); %sacando parches de trabajo
ddHR=patch_ups(ddLR);

%% buscar conjuntos de parches semejantes

B_sim=pch_sim(ddHR,DDHr);

%% Permutaciones

B_per=permutaciones(B_sim);

%% Buscar parches HR con cov semejante
cov_base = extract_cov(cov(ddHR));

[~,~,n_per]=size(B_per);
for k=1:n_per
    covs=extract_cov(cov(B_per(:,:,k)));
    vals(:,k)=dot(cov_base,covs);
end
[m,I]=max(vals);
A=B_per(:,:,I);
A_p=pinv(A);