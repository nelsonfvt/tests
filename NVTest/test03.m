%% Inicio
[rows,cols]=size(DDHr);
t_vecs=find(WNueLi);

%% Upsampling parches LR
ddLR=DDLR(:,t_vecs);
ddHR=zeros(rows,length(t_vecs));

[Xq,Yq,Zq]=meshgrid((1:2/5:3),(1:2/5:3),(1:2/5:3));
for i=1:length(t_vecs)
    pch3d=reshape(ddLR(:,i),[3,3,3]); %parche 3x3x3
    pchint=interp3(pch3d,Xq,Yq,Zq); % parche 6x6x6
    ddHR(:,i)=reshape(pchint,[6*6*6,1]);
    ddHR(:,i)=ddHR(:,i)/norm(ddHR(:,i));
end

%% buscar conjuntos de parches semejantes

dists=zeros(length(t_vecs),cols);
for i=1:length(t_vecs)
    for j=1:cols
        dists(i,j)=dot(ddHR(:,i),DDHr(:,j));
    end
end

%% Buscar parches HR con cov semejante
cov_base = cov(ddHR);