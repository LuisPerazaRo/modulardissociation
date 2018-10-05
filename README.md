# Modular dissociation 
## Data sharing tutorial

Add the tutorial toolbox in the matlab path
```Matlab
addpath E:\Inv_DLB_ModularHypothesis\Datasharing\MDToolbox
```
Add toolbox for the minimum spanning tree
```Matlab
addpath E:\Inv_DLB_ModularHypothesis\MatlabToolboxes\matlab_networks_routines\code
```
Add the the Brain Connectivity Toolbox of Olaf Sporns
```Matlab
addpath E:\Inv_DLB_ModularHypothesis\MatlabToolboxes\BCT\2017_01_15_BCT
```
Load the NKI database matrices
```Matlab
load('NKImatrices_Peraza.mat')
```
Load the coordinates (MNI)
```Matlab
load('MNIcoordinateSystems\MNI100roi_atlas.mat')
coords=MNI100coords(:,1:3);
clear MNI100coords  %delete the variable
```

%Select one participant from the NKI database
pat = 100;
pat_matrix = abs(Connectome100_YA_NKI(:,:,pat)); %The absolute value
deg = 4; % Average node degree (whole network)

%Perform local thresolding for the connectivity matrix
newmatLT=localThresholding(pat_matrix,deg);
newmatLT=newmatLT(:,:,2);
%And estimate its Louvain's modularity
modularityLT=0;
for iter=1:50 
    [auxCiLT,aux_Q] = community_louvain(newmatLT .* pat_matrix);
    if aux_Q > modularityLT
        modularityLT=aux_Q;
        CiLT=auxCiLT;
    end
end

%Perform global thresholding for the connectivity matrix
newmatGT=globalThresholding(pat_matrix,4,'binary');
modularityGT=0;
for iter=1:50 
    [auxCiGT,aux_Q] = community_louvain(newmatGT .* pat_matrix);
    if aux_Q > modularityGT
        modularityGT=aux_Q;
        CiGT=auxCiGT;
    end
end

figure
[X,Y,INDSORT] = grid_communities(CiLT); % call function
imagesc(newmatLT(INDSORT,INDSORT));           % plot ordered adjacency matrix
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'r','linewidth',2);
title('Locally thresholded matrix - communities')
hold off;

figure
[X,Y,INDSORT] = grid_communities(CiGT); % call function
imagesc(newmatGT(INDSORT,INDSORT));           % plot ordered adjacency matrix
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'r','linewidth',2);
title('Globally thresholded matrix - communities')
hold off;

%Compute Modular dissociation - difference between both construction
%methods
MD_LTGT=ModularVariance(CiLT,CiGT);
figure
scatter3(coords(:,1),coords(:,2),coords(:,3),100.^MD_LTGT,'o','fill','MarkerFaceColor','b')
title('Modular Dissociation - Example')


