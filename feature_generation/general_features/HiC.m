clear all;
clc;
load('HiC.mat');
c=readtable('hic_data.csv');
ORF=table2array(c(1:end,1));
mRNA=table2array(c(1:end,2));

thres = 6.58;
d = [];

% remove genes with no data
for i = 1:length(G)
    val_1 = G(i,:);
    if sum(val_1) == 0
        d = [d,i];
    end
end

gene_id(d) = [];
G(d,:) = [];
G(:,d) = [];
ORF2=gene_id;
Mean_dist = zeros(length(G),1);
Median_dist = zeros(length(G),1);
Min_dist = zeros(length(G),1);
rank = zeros(length(G),1);
close_genes_mean_mRNA = zeros(length(G),1);

% create mrna vector according to gene_id
for i=1:length(G)
    z = ORF2{i};
    for j = 1:length(find(contains(ORF,z)))
        o = find(contains(ORF,z));
        if strcmp(ORF{o(j)},z) == 1
            mrna(i,:) = mRNA(o(j),1);
        end
    end
end

% calculate paramters
for i =1:length(G)
    val_1 = G(i,:);
    ind = find(val_1<thres);
    val_2 = val_1(find(val_1<thres));
    ind_0 = find(val_2==0);
    val_2(ind_0) = [];
    ind(ind_0) = [];
    close_genes_mean_mRNA(i)=mean(mrna(ind));
    rank(i) = length(val_2);
    Mean_dist(i)=mean(val_2);
    Median_dist(i)=median(val_2);
    Min_dist(i)=min(val_2);
end

% create the features for all genes according to final ORF, missing values
% filled with nan

for i=1:size(ORF,1)
    check = 0;
    z = ORF{i};
    if sum(contains(ORF2,z))>0
        for j = 1:length(find(contains(ORF2,z)))
            o = find(contains(ORF2,z));
            if strcmp(ORF2{o(j)},z) == 1
              Mean_Dist(i,:) = Mean_dist(o(j),1);
              Median_Dist(i,:) = Median_dist(o(j),1);
              Min_Dist(i,:) = Min_dist(o(j),1);
              close_genes_mean_mrna(i,:) = close_genes_mean_mRNA(o(j),1);
              Rank(i,:) = rank(o(j),1);
              check = 1;
            end
        end
    else
          Mean_Dist(i,1) = nan;
          Median_Dist(i,1) = nan;
          Min_Dist(i,1) = nan;
          close_genes_mean_mrna(i,:) = nan;
          Rank(i,:) = nan;
    end
    if check == 0
          Mean_Dist(i,:) = nan;
          Median_Dist(i,:) = nan;
          Min_Dist(i,1) = nan;
          close_genes_mean_mrna(i,:) = nan;
          Rank(i,:) = nan;
    end
end


final_table = table(ORF,Mean_Dist,Median_Dist,Min_Dist,close_genes_mean_mrna,Rank);
writetable(final_table,'HiC.csv');
