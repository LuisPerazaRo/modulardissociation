# Modular dissociation 
## Data sharing tutorial

Welcome to the data and code sharing wegpage for the article “The functional brain favours segregated modular connectivity at old age unless targeted by neurodegeneration” which is currently under review.

The database for this tutorial is published in [Figshare](https://figshare.com/s/7057a9ac73458c3ebbcc). This database contains a PDF file describing the variables and a Matlab file with all matrices. Download the matlab data file NKImatrices_Peraza.mat and you will be ready for this tutorial!
 
We start adding the Toolboxes we need to the Matlab path. My custom functions from modular variability (MV), and local/global network thresholding are stored in the MDToolbox folder,
```Matlab
addpath MDToolbox
```
Add toolbox for the minimum spanning tree. This functional belongs to the [MTNA Toolbox](http://strategic.mit.edu/downloads.php?page=matlab_networks) and can be downloaded from the link provided
```Matlab
addpath matlab_networks_routines\code
```
Finally, we add the [Brain Connectivity Toolbox, BCT](https://sites.google.com/site/bctnet/) to the Matlab path,
```Matlab
addpath BCT\2017_01_15_BCT
```
Load the NKI database matrices
```Matlab
load('NKImatrices_Peraza.mat')
```
After loading the mat file to the workspace you will see several variables, including base demographics for age and sex and the connectivity matrices. These connectivity matrices were estimated using Pearson correlations and the functional atlas estimated from an independent older adult group from [Peraza et al. 2017](https://onlinelibrary.wiley.com/doi/abs/10.1002/hbm.23499). These atlases are available on [Figshare](https://figshare.com/s/7057a9ac73458c3ebbcc).

Now, download the MNI coordinate system used in my study from the Figshare database and load the ROI coordinates (in MNI). For the example in this tutorial, we will use the 100 ROI atlas and matrices, but the same principles apply to the other atlases.
```Matlab
load('MNIcoordinateSystem\MNI100roi_atlas.mat')
coords=MNI100coords(:,1:3);     %Ignore the fourth columns
clear MNI100coords              %delete the variable and keep coords
```
Select one participant from the NKI database, this participant can be random. It his case we choose 100th participant from the younger adult connectivity matrices, and we specify as well as the average node degre we wish to threshold the connectivity matrix:
```Matlab
pat = 100;
pat_matrix = abs(Connectome100_YA_NKI(:,:,pat)); %The absolute value
deg = 4; % Average node degree (whole network)
```
Remember that in our investigation we took the absolute value of the connectivity matrix. Hence the variable pat_matrix stores the absolute value of the Pearson correlations.

Now we perfom **local thresolding** of this connectivity matrix. Local thresholding is implemented in function localThresholding.mat in the MDToolbox and take the weighted connectivity matrix and the average node degree as parameters. 
```Matlab
newmatLT=localThresholding(pat_matrix,deg);
newmatLT=newmatLT(:,:,2);
```
We also estimate network commnunities with the Louvain's algorithm (from the BCT) after local thresholding. For this example we will run Louvain's 50 times. Current publications recommend > 500 iterations of the algorithm, although this number depends on the networ size. For the code below at each iteration of the for loop, we compare and save the new estimated community stored in CiLT variable that shows a higher modularity index.
```Matlab
modularityLT=0;
for iter=1:50 
    [auxCiLT,aux_Q] = community_louvain(newmatLT .* pat_matrix);
    if aux_Q > modularityLT
        modularityLT=aux_Q;
        CiLT=auxCiLT;
    end
end
``` 
Now we do the same for the **global thresholding** of the connectivity matrix,
```Matlab
newmatGT=globalThresholding(pat_matrix,4,'binary');
modularityGT=0;
for iter=1:50 
    [auxCiGT,aux_Q] = community_louvain(newmatGT .* pat_matrix);
    if aux_Q > modularityGT
        modularityGT=aux_Q;
        CiGT=auxCiGT;
    end
end
```
```Matlab
figure
[X,Y,INDSORT] = grid_communities(CiLT); % call function
imagesc(newmatLT(INDSORT,INDSORT));           % plot ordered adjacency matrix
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'r','linewidth',2);
title('Locally thresholded matrix - communities')
hold off;
```

```Matlab
figure
[X,Y,INDSORT] = grid_communities(CiGT); % call function
imagesc(newmatGT(INDSORT,INDSORT));           % plot ordered adjacency matrix
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'r','linewidth',2);
title('Globally thresholded matrix - communities')
hold off;
```

Compute Modular dissociation - difference between both construction methods
```Matlab
MD_LTGT=ModularVariance(CiLT,CiGT);
figure
scatter3(coords(:,1),coords(:,2),coords(:,3),100.^MD_LTGT,'o','fill','MarkerFaceColor','b')
title('Modular Dissociation - Example')
```

