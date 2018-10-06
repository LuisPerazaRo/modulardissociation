function MV = ModularVariance(Civec1,Civec2)
%This function estimates modular variance between two community definitions
%
%Input:
%Civec1 - Community definition 1
%Civec2 - Community definition 2
%
%Output:
%MV - Vector with modular variability estimates
%
%It is necesary for this function that Civec1 and Civec2 are of the same
%length. It is not necesary that both community definitions have the same 
%number of modules. 
%
%This function is part of the research by Peraza et al.
%2018 "The functional brain favours segregated modular connectivity at old
%age unless targeted by neurodegeneration".

Civec1=Civec1(:);
Civec2=Civec2(:);

if numel(Civec1)~=numel(Civec2)
   error('Modularity definitions are not of same length')
end

MV=zeros(size(Civec1));
for iter=1:numel(Civec1)
    mod1_ind=Civec1==Civec1(iter);
    mod2_ind=Civec2==Civec2(iter);
    
    modunion=sum(mod1_ind.*mod2_ind);
    mod1_nodes=sum(mod1_ind);
    mod2_nodes=sum(mod2_ind);
    
    MV(iter)=1-(modunion^2)/(mod1_nodes*mod2_nodes);
end