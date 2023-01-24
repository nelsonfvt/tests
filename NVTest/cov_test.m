% Experimento covarianzas de diccionarios
% Cargar datos

load('observ.mat')
load('weights.mat')
load('diccor.mat')
load('diccmd.mat')

%Dr_HR random dictionaries
%Dr_HRn random dictionaries + noise
%Dx_HR upsampled dictionaries

%HR_arr parches HR originales
%LR_arr parches LR originales

%HRnear_arr 
%LRnear_arr

%W pesos para LR
%Whr pesos para HR

%% selecionar un atomo

atm = 890;
%parche
HR_ptch = HR_arr(:,atm);
LR_ptch = LR_arr(:,atm);
%pesos
Wp = W(:,atm);
Whrp = Whr(:,atm);
%extrae atomos diccs
idx = find(Wp);
Wp = Wp(idx);
Dnear_lr = LRnear_arr(:,idx);
Dnear_hr = HRnear_arr(:,idx);
Dx_HRo = Dx_HR(:,idx,:);
Dr_HRo = Dr_HR(:,idx,:);
Dr_HRno = Dr_HRn(:,idx,:);
%atomos dicc HRhr
idx = find(Whrp);
Whrp = Whrp(idx);
D_hro = HRnear_arr(:,idx);

%organizar pesos y atomos
[Wp,idx] = sort(abs(Wp),'descend');
Dnear_lr = Dnear_lr(:,idx);
Dnear_hr = Dnear_hr(:,idx);
Dx_HRo = Dx_HRo(:,idx,:);
Dr_HRo = Dr_HRo(:,idx,:);
Dr_HRno = Dr_HRno(:,idx,:);
[Whrp,idx] = sort(abs(Whrp),'descend');
D_hro = D_hro(:,idx);

% extraer filas de diccs
fp = 14;
fps = [87 88 93 94 123 124 129 130];
fila_lr = Dnear_lr(fp,:);
fila_hr = Dnear_hr(fps,:);
filax_HRo = Dx_HRo(fps,:,:);
filar_HRo = Dr_HRo(fps,:,:);
filar_HRno = Dr_HRno(fps,:,:);
fila_hro = D_hro(fps,:);

% covs

cov_hr = Covs(fila_lr,fila_hr);
cov_HRo = Covs(fila_lr,filax_HRo);
covr_HRo = Covs(fila_lr, filar_HRo);
covr_HRno = Covs(fila_lr, filar_HRno);
cov_hro = Covs(fila_lr, fila_hro);
%
tt = Covs(fila_lr,fila_lr);
aaa = covr_HRo - covr_HRno;
mean(mean(aaa))

%% funciónes auxiliares

function [Mx] = Covs(flr,fsx)
    tm = size(fsx);
    if length(tm) > 2
        nf = tm(1);
        nc = tm(3);
        Mx = zeros(nf,nc);
        for i=1:nc
            for j=1:nf
                Mx(j,i) = trace(mx_covar(flr,fsx(j,:,i)));
            end
        end
    else
        Mx = zeros(tm(1),1);
        for i=1:tm(1)
            Mx(i) = trace(mx_covar(flr,fsx(i,:)));
        end
    end
end

function [val] = medida (atm, W, DLR, DHR)
    % seleciona pesos para parches atm
    Wp = W(:,atm);
    
    %extrae atomos diccs
    idx = find(Wp);
    Wp = Wp(idx);
    DLR_s = DLR(:,idx,:);
    DHR_s = DHR(:,idx,:);
    
    %organizar pesos y atomos
    [Wp,idx] = sort(abs(Wp),'descend');



    % extraer filas de diccs
    fp = 14;
    fps = [87 88 93 94 123 124 129 130];

end