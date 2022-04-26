%% Inicio
[rows,cols]=size(DDHr);
t_vecs=find(WNueLi);

%% Upsampling parches LR
ddLR=DDLR(:,t_vecs);
ddHR=zeros(rows,length(t_vecs));

[Xq,Yq,Zq]=meshgrid((1:2/5:3),(1:2/5:3),(1:2/5:3));
for i=1:length(t_vecs)
    pch3d=reshape(ddLR(:,i),[3,3,3]); %parche 3x3x3
    pchint=interp3(pch3d,Xq,Yq,Zq); %interpolando a parche 6x6x6
    ddHR(:,i)=reshape(pchint,[6*6*6,1]);
    ddHR(:,i)=ddHR(:,i)/norm(ddHR(:,i)); %vectores unitarios
end

%% buscar conjuntos de parches semejantes

n_sim=5;
dists=zeros(length(t_vecs),cols);
B_sim=zeros(rows,length(t_vecs),n_sim);
for i=1:length(t_vecs)
    for j=1:cols
        dists(i,j)=dot(ddHR(:,i),DDHr(:,j));
    end
    [B,I]=sort(dists(i,:),'descend');% ordenando desde el más parecido
    for k=1:n_sim
        B_sim(:,i,k)=DDHr(:,I(k));% llenando arreglo con más parecidos
    end
end

%% Buscar parches HR con cov semejante
cov_base = extract_cov(cov(ddHR));

for k=1:n_sim
    covs(:,k)=extract_cov(cov(B_sim(:,:,k)));
    vals(:,k)=dot(cov_base,covs(:,k));
end
