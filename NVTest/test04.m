%% Inicio
parches=Dhr10*WHr10; %parches originales HR
%Dlr10=patch_dws(Dhr10); % generando dicc LR
[p_rows,p_cols]=size(parches);
[rows,cols]=size(Dhr10);
t_vecs=find(WHr10(:,1));
n_atm=length(t_vecs);

%% Upsampling parches LR
ddHr = patch_ups(Dlr10);

%% Encontrar el que tiene la covarianza mas parecida
% Uno por uno de los parches del volumen
A=zeros(rows,n_atm,p_cols);
I=zeros(1,p_cols);

for p=1:p_cols %numero de parches del volumen
    p_idx=find(WHr10(:,p));% posiciones de parches
    pcHR=ddHr(:,p_idx); %atomos sobremuestreados
    B_sim=pch_sim(pcHR,Dhr10); %buscando atomos similares
    B_per=permutaciones(B_sim); %permutaciones
    vals = covars(pcHR,B_per); %covarianzas
    [~,I(p)]=min(vals); %min porque son dist euclid.
    A(:,:,p)=B_per(:,:,I(p));
end

%% Reconstrucciones

rpatchs=zeros(size(parches));
for p=1:p_cols
    ws=WHr10(:,p);
    ws_id=find(ws);
    ws=ws(ws_id);
    rpatchs(:,p)=A(:,:,p)*ws;
end