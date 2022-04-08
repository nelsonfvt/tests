%% Inicio
[rows,cols]=size(DDHr);
t_vecs=find(WNueLi);

%% Submuestreo media local (paper)
ddLR=DDLR;

%%Test
% extract vectors
ddLRs=ddLR(:,t_vecs);
ddHRs=DDHr(:,t_vecs);
% normalize LR atoms
for i=1:length(t_vecs)
    ddLRs(:,i)=ddLRs(:,i)/norm(ddLRs(:,i));
end

%uno contra todos
clr=[];
ul=ddLRs(:,1);
%covarianzas de atomos LR
for i=2:length(t_vecs)
    vl=DDLR(:,i);
    aa=[ul vl];
    c=cov(aa);
    clr(i-1)=c(2,1);
end
%covarianzas entre todos los atomos HR
uh=ddHRs(:,1);
chr=[];
for i=2:cols
    vh=DDHr(:,i);
    aa=[uh vh];
    c=cov(aa);
    chr(i-1)=c(2,1);
end

%Encontrar las covarianzas mas parecidas
nf=length(clr);
nc=length(chr);
Arr=zeros(nf,nc);
for f=1:nf
    for c=1:nc
        Arr(f,c)=abs(clr(f)-chr(c));
    end
end
[M,I]=min(Arr,[],2);
