%% Inicio
[rows,cols]=size(DDHr);
t_vecs=find(WNueLi);

%% Submuestreo impares
ps=1:2:rows;
ddLR=DDHr(ps,:);

%% Submuestreo pares
ps=2:2:rows+1;
ddLR=DDHr(ps,:);

%% Submuestreo de 8
ps=1:8:rows;
ddLR=DDHr(ps,:);

%% Submuestreo de 8
ps=3:8:rows;
ddLR=DDHr(ps,:);

%% Submuestreo media local (paper)
ddLR=DDLR;

%% Test
% extract vectors
ddLRs=ddLR(:,t_vecs);
ddHRs=DDHr(:,t_vecs);
% normalize LR atoms
for i=1:length(t_vecs)
    ddLRs(:,i)=ddLRs(:,i)/norm(ddLRs(:,i));
end

% 1 contra todos
res=[];
ul=ddLRs(1,:);
uh=ddHRs(1,:);
for i=2:length(t_vecs)
    %LR
    vl=ddLRs(i,:);
    ct=dot(ul,vl)/(norm(ul)*norm(vl));
    res(1,i-1)=acosd(ct);
    %HR
    vh=ddHRs(i,:);
    ct=dot(uh,vh)/(norm(uh)*norm(vh));
    res(2,i-1)=acosd(ct);
end