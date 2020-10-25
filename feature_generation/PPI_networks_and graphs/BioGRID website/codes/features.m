clc
clear all
T = readtable('try.csv', 'HeaderLines',1);  
Systematic_A=T{1:end,1};
Systematic_B=T{1:end,2};
out1 = unique(Systematic_A,'stable');
out2 = unique(Systematic_B,'stable');
all_proteins=[out1;out2];
all_proteins_unique=unique(all_proteins,'stable');
rank=zeros(1,length(all_proteins_unique));
for i=1:length(out1)
    for j=1:length(Systematic_B)
        if isequal(Systematic_B{j},out1{i})
            Systematic_B{j}=Systematic_A{j};
            Systematic_A{j}=out1{i};
        end
    end
   rank(j)=sum(strcmp(out1{i},unique(Systematic_A)));
end 
% app2 = strcmp('YOR123C',Systematic_B);
% 
%     

