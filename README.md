# Modular dissociation 
## Data sharing tutorial

Welcome to the data and code sharing tutorial for the article “The functional brain favours segregated modular connectivity at old age unless targeted by neurodegeneration” which is currently under review.

The database for this tutorial is published in [Figshare](https://figshare.com/s/7057a9ac73458c3ebbcc). This database contains a PDF file describing the variables as well as a Matlab file with all participant matrices from the NKI database. Download the matlab data file NKImatrices_Peraza.mat and you will be ready for this tutorial!
 
We start adding the toolboxes we need to the Matlab path. My custom functions for modular variability/dissociation (MV/MD), and local/global network thresholding are stored in the MDToolbox folder,
```Matlab
addpath MDToolbox
```
Also, add the toolbox for the minimum spanning tree (MST, needed for local thresholding). The function that estimates the MST belongs to the [MTNA Toolbox](http://strategic.mit.edu/downloads.php?page=matlab_networks) and can be downloaded from the link provided,
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
After loading the mat file to the Matlab workspace you will see several variables. These comprise basic demographics for age and sex, and the connectivity matrices. These connectivity matrices were estimated using Pearson correlations and the functional atlases from an independent older adult group from [Peraza et al. 2017](https://onlinelibrary.wiley.com/doi/abs/10.1002/hbm.23499). These atlases are available on the [Figshare](https://figshare.com/s/7057a9ac73458c3ebbcc) repository.

Now, download the MNI coordinate system used in my study from the Figshare database and load the ROI coordinates (in MNI). For the example in this tutorial, we will use the 100 ROI atlas and matrices, but the same principles apply to the other atlases.
```Matlab
load('MNIcoordinateSystem\MNI100roi_atlas.mat')
coords=MNI100coords(:,1:3);     %Ignore the fourth columns
clear MNI100coords              %delete the variable and keep coords
```
Select one participant from the NKI database, and this participant can be random for the purposes of this tutorial. It his case I choose the 100th participant from the younger adult connectivity matrices, and we specify as well as the average node degre we wish to threshold the connectivity matrix:
```Matlab
pat = 100;
pat_matrix = abs(Connectome100_YA_NKI(:,:,pat)); %The absolute value
deg = 4;                                         % Average node degree (whole network)
```
Remember that in this investigation I took the absolute value of the connectivity matrix. Hence, the variable pat_matrix stores the absolute value of the Pearson correlations.

Now we perfom **local thresolding** of this connectivity matrix. Local thresholding is implemented in function localThresholding.mat within the MDToolbox, and takes the weighted connectivity matrix and the average node degree as parameters. 
```Matlab
newmatLT=localThresholding(pat_matrix,deg);
newmatLT=newmatLT(:,:,2);
```
With the thresholded matrix, estimate network commnunities with the Louvain's algorithm (from the BCT). For this example I will run Louvain's 50 times. Current publications recommend > 500 iterations of the algorithm, although this number depends on the network size. For the code below, I am comparing and saving, at each iteration of the for loop, the new estimated community stored in CiLT variable that shows a higher modularity index Q.
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
Now let's do the same for the **global thresholding** of the connectivity matrix, and save the estimated community definition in the CiGT vector.
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
Now, inspect the estimated community defnitions with the following code and by borrowing some functions from the BCT toolbox. Plot the community matrices for the local and global threshold approaches.

```Matlab
figure
[X,Y,INDSORT] = grid_communities(CiLT);   % call function
imagesc(newmatLT(INDSORT,INDSORT));       % plot ordered adjacency matrix
hold on;                                  % hold on to overlay community visualization
plot(X,Y,'r','linewidth',2);
title('Locally thresholded matrix - communities')
hold off;
```
![Local Threshold](/images/LocalThreshold_Example.png)

```Matlab
figure
[X,Y,INDSORT] = grid_communities(CiGT); % call function
imagesc(newmatGT(INDSORT,INDSORT));           % plot ordered adjacency matrix
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'r','linewidth',2);
title('Globally thresholded matrix - communities')
hold off;
```
![Global Threshold](/images/GlobalThreshold_Example.png)

Notice how local and global thresholding methods led to different community definitions. When the connectivity matrix is thresholded **locally** communities are less sparse and larger, while in **global threshold** few communities concentrate the mayority of the strongest links. This is because local thresholding is based on the _k_ nearest neighbour graph (_k_-NNG) which favours creation and segregation of communities. Remember as well that modularity Q is larger when matrices are localy thresholded.

With these two definitions, we can compute the modular dissociation/variability statistic, which is defined by

![MD definition](/images/MD_definition.png)

where the numerators represent the number of nodes in common between communities X(i) and X(j), and the denominators represent the number of nodes within the community. Here node s belongs to the network i and similar for network j.

This definition of Modular Dissociation (MD) is implemented in the MDToolbox function ModularVariance.mat and it is executed as follows:

```Matlab
MD_LTGT=ModularVariance(CiLT,CiGT);
```
We can now see the nodal MD values with a scatter3 plot in Matlab, and using the MNI coordinate definition for the 100-ROI atlas,

```Matlab
figure
scatter3(coords(:,1),coords(:,2),coords(:,3),100.^MD_LTGT,'o','fill','MarkerFaceColor','b')
title('Modular Dissociation - Example')
```
![MD example](/images/MD_example.png)

Notice the brain shape of the nodal arrangement, and how the motor-sensory, occipital and temporal cortices show, for this participant, low values of MD while high values are present in the frontal cortex.

I hope you liked this tutorial and don't forget to visit our lab webpage www.lewybodylab.org, which has further information about Lewy body dementia research.

Best.
