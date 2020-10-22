
clear all
close all
clc
%cd('C:\Users\karinsio\Documents\המסמכים שלי\לימודים\שנה ד\סמסטר ב\iGEM\מודל חדש\Initial Feature Selection CFS') %%% change datapath.
%% 
%%% Update 16.08: 
%%% 1. When using "max" on the correlation matrices, to find the first
%%% feature, I added "abs". It didnt change the outcome but still needed to be fixed.  
%%% 2. Only GFP is considered. The RFP is for "transfer learning". 

%%% Update 22.09:
%%% 1. Target Gene features (longest common sub string, profile features). 
%%% 2. Conservation features (genes with missing values where mostly filtered - about 4700 genes left). 

%% load data

% extract features:
features_tbl = readtable('normalized_updated_new_features_table.csv');
% % add GFP cARS: 
% GFP_cARS = readtable('GFP_CARS.csv');
% [~,unique_ind] = unique(GFP_cARS(:,1));
% GFP_cARS = GFP_cARS(unique_ind,:); %unique indices
% joined_table = join(features_tbl,GFP_cARS); %join the features tables
% cars = table2array(joined_table(:,end)); %for normalization
% joined_table(:,end)=num2cell((cars-mean(cars))/std(cars)); %normalization
% features_tbl=joined_table; %back to original name

% unique values only:
[~,unique_ind] = unique(features_tbl(:,1));
X = table2array(features_tbl(unique_ind,2:end)); % feature matrix, without 'ORF' column. 

% extract labels for GFP only (the RFP is for "transfer learning"):
labels_folder = [pwd filesep 'flourescence_tables (September)'];
tables_names = ["flourescence_table_NATIVEpr_GFP.csv", "flourescence_table_NOP1pr_GFP.csv"];

% y = zeros(size(X,1),length(tables_names)); 
y = [];
for i=1:length(tables_names)
    label_tbl = readtable([labels_folder filesep char(tables_names(i))]);
    y = [y, table2array(label_tbl(:,2:end))]; 
end

% fill nan values in the labels columns:
for j=1:size(y,2)
    y(isnan(y(:,j)),j) = mean(y(~isnan(y(:,j)),j)); %fill nans with mean value. 
end

feature_feature_correlation = corr(X,'type','Spearman'); 
feature_class_correlation = corr(X, y, 'type', 'Spearman'); 
[~,max_ind] = max(abs(feature_class_correlation)); %output size is [1,size(y,2)].
max_ind = mode(max_ind); %the first chosen feature is the most correlative with all the labels 

%% Feature selection: SFS with CFS (filter method)
N = 300; %the number of chosen features
for n_features=2:N
    %Looking for the best n_features combination
    CFS=zeros(size(X-length(max_ind),2),1);
    available_features=setdiff((1:size(X,2)),max_ind); %all the features beside the chosen ones
    for r=available_features
        CFS(r)=calculate_CFS(feature_class_correlation,feature_feature_correlation,[max_ind,r]);
        %disp(['CFS criterion of feature(s) ',num2str(max_ind),' with feature #',num2str(r),' = ',num2str(CFS(r))])
    end
    [~,ind]=max(CFS); %chosen feature in this step
    max_ind = [max_ind,mode(ind)]; %add to the list of chosen features indices
end
chosen_features = max_ind; %indices of chosen features (columns in X matrix).

%% chosen features names: 
feature_names = features_tbl.Properties.VariableNames; %column names of 'features_tbl'
feature_names = feature_names(2:end); %without 'ORF' because it is not a feature
chosen_features_names = feature_names(chosen_features)';

%% save file:
% the file will contain the feature names and indices (matlab-wise). 
output_tbl = [chosen_features_names, num2cell(chosen_features)'];
output_tbl = cell2table(output_tbl,'VariableNames',{'feature_name','feature_index_matlab'});
writetable(output_tbl,'chosen_features.xlsx');
