function mat = globalThresholding(Cmat,deg,Mtype)
%Threshold the undirected connectivity matrix given a desired 
%average node degree.
%
%Usage: 
%   mat = globalThresholding(Cmat,deg,Mtype)
%
%Output: 
%   mat - weighted directed or binary matrix
%
%Input:
%   Cmat - weighted symetric/undirected connectivity matrix
%   deg - desired average network node degree in mat
%   Mtype - string variable: Mtype = 'weighted','binary'

[re,co]=size(Cmat);
if re~=co
    error('matrix is not square');
end

%Use the upper triangular absolute part only
Cmat=abs(triu(Cmat,1));

edge1=round(deg*re/2);
vec=reshape(Cmat,1,re*co);
[vecsort,ind]=sort(vec,2,'descend');
vecsort(edge1+1:end)=0;
if strcmp(Mtype,'binary')==1
    vecsort(1:edge1)=1;
elseif strcmp(Mtype,'weighted')==1
    vecsort(1:edge1); %do nothing
else
    disp('Warning: undefined matrix type')
end

vec(ind)=vecsort;
mat=reshape(vec,re,co);

%Make the matrix symmetric again
mat=mat+mat';