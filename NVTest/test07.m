%Experimento para demostrar Proc. Gauss. Estacionario

%Se requiere usar Nifti_tools
%addpath('/usr/local/MATLAB/matlab_comm/NIfTI_tools/')
%addpath('/home/nelsonvt/Data_Files/Harmonization/GE_Atest/Cases/Matlab_DWI/')

%%Cargando datos DWI
gen_path = ['/home/nelsonvt/Data_Files/Harmonization/GE_Atest/Cases/A/GE/'];
stats = my_preprocess(gen_path);


%% Calculo diferencias
ix=find(stats.cluster==7); % Seleccionando un cluster
%stats.vox(ix,:) vectores intensidades x cda voxel
%dc_arr = my_diffcovar(stats.vox(ix,:));

%calcular diferencias (vectores)
difs = my_difers(stats.vox(ix,:));
clear ix

%% Generar diferentes configuraciones (voxeles)
% Se toman dos voxeles diferentes (aleatorios)
clc
i_vx = randperm(length(difs(:,1,1)),2); 
%[8178,8625]; SI
%[7127,71]; NO
vx1 = reshape(difs(i_vx(1),:,:), [30 30]);
vx2 = reshape(difs(i_vx(2),:,:), [30 30]);

%% Generar diferentes configuraciones (ex gradientes)
% sacar los mismos gradientes aleatorios (filas y columnas)
i_gr = sort(randperm(30,1));
D1 = vx1;
D1(i_gr,:) = [];
D1(:,i_gr) = [];

D2 = vx2;
D2(i_gr,:) = [];
D2(:,i_gr) = [];

%% Calculo covarianzas y comparación
%C = my_covars(difs);
D1 = D1(1,:);
D2 = D2(1,:);
C1 = my_covar(D1);
C2 = my_covar(D2);
 
% generando descomposición matricial
[X1,s1,~] = svd(C1);
[X2,s2,~] = svd(C2);

% Comparación Matrices de covarianza
V11 = D1*X1(:,1);
V22 = D2*X2(:,1);
V12 = D1*X2(:,1);
V21 = D2*X1(:,1);

% cálculo S-statistics
S1 = 2*sum((V11-V21).^2 + (V12-V22).^2);
S2 = sum(((V11+V22) - (V12+V21)).^2);
S3 = sum(((V11+V12) - (V21+V22)).^2);

% prod punto Vector principal
d_pr = X1(:,1)' * X2(:,1);

formatSpec = 'Valor %7.4f\n';
fprintf(['S1: ' formatSpec],S1)
fprintf(['S2: ' formatSpec],S2)
fprintf(['S3: ' formatSpec],S3)
fprintf(['Pr_p: ' formatSpec],d_pr)

%% Funciones auxiliares
function [stats] = my_preprocess(gen_path)
%Funcion para realizar la carga de datos y preprocesamiento:
%   
%Input:
%   path: cadena con ruta a los datos
%Output:
%   stats: estructura con los datos agrupados

%Cargando imagenes DWI
dwi=load_untouch_nii([gen_path, 'bet_st/dwi_bet.nii']);
bval=load([gen_path,'bet_st/dwi.bval']);
bvecs=load([gen_path, 'bet_st/dwi.bvec']);

%excluyendo b0s
ps=find(bval);
bst=dwi.img(:,:,:,ps); % gradientes
bvec=bvecs(ps,:); % direcciones

Nclust=15;
stats=ex_statvoxels(1,bst);

% Clusters
TM_model.C=[];
TM_model.range=[];
TM_model.ws=[];
TM_model.fits={};
Xo=[stats.voxmean stats.voxstd];
%kmeans y graficas
[GMidx,gmC]=kmeans(Xo,Nclust,'MaxIter',1000,'Options',statset('UseParallel',1),'Start','cluster','Replicates',5);

TM_model.C=gmC;
%ordenar datos
[~,indm]=sortrows(TM_model.C);

TM_model.C=TM_model.C(indm,:); %ordenando centros
stats.cluster=zeros(size(GMidx));

for c=1:Nclust
    nc_idx=find(indm==c); % Encontrando el nuevo indice del cluster c
    n_pos=find(GMidx==c); % Encontrar donde están los puntos correspondeintes al cluster
    stats.cluster(n_pos)=nc_idx;
end

end

function [difs]=my_difers(vectores)
%Función para calcular todas las diferencias de intensidades
%Input:
%   vectores: Matriz con vectores de intensidades
%Output:
%   difs: arreglo 3D con diferencias
[fil,col] = size(vectores);
difs = zeros(fil,col,col);

for k = 1:fil %
    for i = 1:col % 
        for j = 1:col
            difs(k,i,j) = vectores(k,i)-vectores(k,j);
        end
    end
end
end

function [Mx] = my_covar(Arr)
%Función para calcular la matriz covarianza de un vector
%Input:
%   Arr: Vector fila [1, N]
%Output:
%   Mx: Matriz de covarianza
m=mean(Arr);
At=Arr'-m;
A=Arr-m;
Mx=(At*A)./length(A)-1;
end
