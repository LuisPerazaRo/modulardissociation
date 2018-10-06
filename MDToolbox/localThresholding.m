function MST_Kmatrix = localThresholding(auxmat,desiG_vec)
%
%This function threshold matrices using local thresholding method.
%Receives a weighted connectivity matrix and a vector with the desired 
%average node degrees at which the user wants to threshold the matrix.
%
%Returns binary undirected matrices in a 3D matrix - the 3rd dimension
%corresponds to the thresholds specified in the second input parameter.
%
%Usage:
%
% MST_Kmatrix = localThresholding(WU_matrix,degree_vec);
% Where:
% WUmatrix - is a wighted undirected matrix
% degree_vec = Is a vector with desired degrees, eg: vec = 3:10
%
%The returned 3D matrix will have [nodes x nodes x length(vec)+1]
%dimensions where +1 is the original minimum spanning tree (MST).

%First elimiante the diagnonal in the weithed matrix
nodes=size(auxmat,2); %Total number of nodes in the matrices
auxmat=auxmat.*(1-diag(ones(nodes,1)));

%Save all the N matrices, always save the primary MST 
MST_Kmatrix=zeros(nodes,nodes, length(desiG_vec)); 

%Prepare the seed MST network from the connectivity matrix. This part uses
%the min_span_tree function from the MIT Matlab Tools for Network analysis 
%toolbox (http://strategic.mit.edu/downloads.php?page=matlab_networks).

Kmatrix=min_span_tree(1-auxmat);
submat2=auxmat; %The Nearest Neighbour matrix to be ubdated

%Use this network as seed for the algorithm
submat2(Kmatrix==1)=0;
%The initial network degree given by the MST
deg=sum(sum(Kmatrix))/(nodes); %The current Kmatrix average degree

if deg > desiG_vec(1) %Error if the desired start is lower than the minimun posible
   error('Desired minimum degree is lower than the MST one')
end

count=0; % Count for security break
%Save the starting MST network
MST_Kmatrix(:,:,1)=Kmatrix;
mindex=1; %Start the matrix index for saving matrices
while deg<=desiG_vec(end) 
    [vals,ind]=max(submat2,[],1);
    [~,sind]=sort(vals,'descend');
    
    if mindex > length(desiG_vec) %Do not break the program because of mindex
        mindex=length(desiG_vec);
    end
    
    for iter=sind
        if sum(Kmatrix(iter,:))<desiG_vec(mindex);
            Kmatrix(iter,ind(iter))=1;
            Kmatrix(ind(iter),iter)=1;
            submat2(iter,ind(iter))=0;
            submat2(ind(iter),iter)=0;       
        end
        
        deg=sum(sum(Kmatrix))/(nodes);
        MST_Kmatrix(:,:,mindex+1)=Kmatrix;
         
        if deg>=desiG_vec(end);
           disp('Degree achieved')
           break; %Break the for 
        end
    end
   
   %security break 
   count=count+1;
    if count==10000 || deg>=desiG_vec(end);
        disp('Break activated')   
        break;
    end
end
