clear all;
clc;
load('CoExp.mat');
A={};
M=[];
max_score=max(max(G));
Weighted_Rank=[];
for i=1:6664
    A{i,1}=G(i,(~isnan(G(i,:))));
    M(i,1)=length(A{i,1});
    Weighted_Rank(i,1)=abs(sum(A{i,1})/(max_score*M(i,1)));
end
ORF=gene_id;
final_table=table(ORF,Weighted_Rank);
%writetable(final_table,'Weighted_Rank.csv');
