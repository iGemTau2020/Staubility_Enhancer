% searching ranks of proteins with updated graph

f = fopen('4932.protein.links.v11.0.txt');
PPI = textscan(f,'%s %s %f','HeaderLines',1); %EndNodes_1, EndNodes_2, Weight
fclose(f);

%extracting proteins names and weigth of interaction
interactions = cellfun(@(x) x(6:end), [PPI{1,1:2}],'UniformOutput',false);% we want the ORF which begins after "4932." i.e from the 6'th char
weight = [PPI{1,3}];
max_weight = max(weight);
ORF_proteins = unique(interactions(:,1));
clear PPI

ranks = zeros(size(ORF_proteins));
weighted_rank = zeros(size(ORF_proteins));
lines_num = length(ORF_proteins);

for i=1:size(ORF_proteins,1)
    if mod(i,10) == 0
        fprintf("at %d out of %d\n", i, lines_num);
    end
    protein = ORF_proteins{i};
    indices_logical = contains(interactions(:,1),protein); %logical vector with same length as interctions vector
    rank_i = length(find(indices_logical)); %how many interactions does the protein have
    ranks(i) = rank_i;
    weighted_rank(i) = sum(double(indices_logical).* weight)/ (max_weight * rank_i);% mean of the weigth for each gene
    %divided by the value of the max weight in the graph
    
end

save('ranks.mat','ranks');
save('ORF_proteins.mat','ORF_proteins');
save('weighted_rank.mat','weighted_rank')
tbl = cell2table(ORF_proteins);
tbl2 = array2table(ranks);
tbl3 = array2table(weighted_rank);
tbl_csv = [tbl,tbl2,tbl3];
writetable(tbl_csv,'rank_data_updated.csv'); 