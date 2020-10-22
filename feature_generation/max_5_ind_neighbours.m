clc;
f = fopen('4932.protein.links.v11.0.txt');
PPI = textscan(f,'%s %s %f','HeaderLines',1); %EndNodes_1, EndNodes_2, Weight
fclose(f);

%extracting proteins names and weigth of interaction
interactions = cellfun(@(x) x(6:end), [PPI{1,1:2}],'UniformOutput',false);% we want the ORF which begins after "4932." i.e from the 6'th char
weight = [PPI{1,3}];
max_weight = max(weight);
ORF_proteins = unique(interactions(:,1));
clear PPI

% mat of ORF index
c = readtable('seq_orf.csv');
ORF = table2array(c(1:end,1));
indexes = zeros(1,6574);
indexes = indexes -1;

for i=1:size(ORF,1)
    z=ORF{i};
    if sum(contains(ORF_proteins,z))>0
        for j=1:length(find(contains(ORF_proteins,z)))
            o=find(contains(ORF_proteins,z));
            if strcmp(ORF_proteins{o(j)},z)==1
                indexes(o(j)) = i;
            end
        end
    end 
end


max_5_ind = zeros(6711,5);
max_5_ind = max_5_ind-1;


for i=1:size(ORF_proteins,1)
    if i~=2000 && i ~=1999
        real_ind = indexes(i);
        protein = ORF_proteins{i};
        indices_logical = contains(interactions(:,1),protein); %logical vector with same length as interctions vector
        weights = weight(find(indices_logical));
        proteins = interactions(:,2);
        proteins = proteins(find(indices_logical));
        [weights_sorted,k] = sort(weights,'descend');
        x = length(weights_sorted);
        if ~isempty(weights_sorted)
            for j = 1:min(5,x)
                p = proteins(k(j));   
                o=find(contains(ORF_proteins,p));
                for h=1:length(o)
                    if strcmp(ORF_proteins{o(h)},p)==1
                        m_i = indexes(o(h));
                    end 
                end
                max_5_ind(real_ind,j)=m_i;
            end
        end
    end
end

max_5_ind = max_5_ind -1;
writetable(table(max_5_ind),'max_5_ind.csv');


