clc
clear all
T = readtable('biogrid-yeastyeast.csv', 'HeaderLines',1);  
Systematic_A=T{1:end,1};
Systematic_B=T{1:end,2};
Experimental_Type=T{1:end,6};
Experimental_Name=T{1:end,5};
Interaction_Throughput=T{1:end,9};
unique_interactions=containers.Map('KeyType','char','ValueType','double');%creating a dictionary
for i=1:length(Systematic_A)
    key=[string(Systematic_B{i}),string(Systematic_A{i}),string(Experimental_Name{i})];
    key=char(join(key,"-"));
    key2=[string(Systematic_A{i}),string(Systematic_B{i}),string(Experimental_Name{i})];
    key2=char(join(key2,"-"));
    if isKey(unique_interactions,key) %if B-A with same experimental type is already a key, delete it
        Systematic_A{i}=[];
        Systematic_B{i}=[];
        Experimental_Type{i}=[];
        Experimental_Name{i}=[];
        Interaction_Throughput{i}=[];
    elseif isKey(unique_interactions,key2) %if A-B with same experimental type is alreay a key, delete it
        Systematic_A{i}=[];
        Systematic_B{i}=[];
        Experimental_Type{i}=[];
        Experimental_Name{i}=[];
        Interaction_Throughput{i}=[];
    else %otherwise, the interaction is new and count it
        key=[string(Systematic_A{i}),string(Systematic_B{i}),string(Experimental_Type{i})];
        key=char(join(key,"-"));
        unique_interactions(key)=1;
    end 
end 
%remove the empty cells from the columns
Systematic_A=Systematic_A(~cellfun('isempty',Systematic_A));
Systematic_B=Systematic_B(~cellfun('isempty',Systematic_B));
Experimental_Type= Experimental_Type(~cellfun('isempty', Experimental_Type));
Experimental_Name=Experimental_Name(~cellfun('isempty',Experimental_Name));
Interaction_Throughput=Interaction_Throughput(~cellfun('isempty',Interaction_Throughput));

%creating a united genes list without duplications
Systematic_A_uni=unique(Systematic_A,'stable');
Systematic_B_uni=unique(Systematic_B,'stable');
all_proteins=[Systematic_A_uni;Systematic_B_uni];
all_proteins_unique=unique(all_proteins,'stable');

% %initalzing all the peatures we would like to calculate
% rank=zeros(1,length(all_proteins_unique));
% num_of_genetic_interactions=zeros(1,length(all_proteins_unique));
% num_of_physical_interactions=zeros(1,length(all_proteins_unique));
% exp_names=unique(Experimental_Name,'stable');
% type_of_interactions=zeros(length(all_proteins_unique),length(exp_names));
% num_of_high=zeros(1,length(all_proteins_unique));
% num_of_low=zeros(1,length(all_proteins_unique));
% num_of_both=zeros(1,length(all_proteins_unique));
% %calculating the features
% for i=1:length(all_proteins_unique)
%        rank(i)=sum(strcmp(all_proteins_unique{i}, Systematic_A))+sum(strcmp(all_proteins_unique{i}, Systematic_B));
%        index_inA=find(strcmp(all_proteins_unique{i}, Systematic_A));
%        index_inB=find(strcmp(all_proteins_unique{i}, Systematic_B));
%        for j=1:length(index_inA)
%            if strcmp(Experimental_Type{index_inA(j)},"genetic")
%                num_of_genetic_interactions(i)=num_of_genetic_interactions(i)+1;
%            else
%                num_of_physical_interactions(i)=num_of_physical_interactions(i)+1;
%  
%            end 
%        end 
%        for j=1:length(index_inB)
%            if strcmp(Experimental_Type{index_inB(j)},"genetic")
%                num_of_genetic_interactions(i)=num_of_genetic_interactions(i)+1;
%            else
%                num_of_physical_interactions(i)=num_of_physical_interactions(i)+1;
%  
%            end
%        end
%        for j=1:length(index_inA)
%           interaction_type=find(strcmp(Experimental_Name{index_inA(j)},exp_names));
%           type_of_interactions(i,interaction_type)= type_of_interactions(i,interaction_type)+1;
%        end
%        for j=1:length(index_inB)
%           interaction_type=find(strcmp(Experimental_Name{index_inB(j)},exp_names));
%           type_of_interactions(i,interaction_type)= type_of_interactions(i,interaction_type)+1;
%        end 
%        for j=1:length(index_inA)
%            if strcmp(Interaction_Throughput{index_inA(j)},"High Throughput")
%                num_of_high(i)=num_of_high(i)+1;
%            elseif strcmp(Interaction_Throughput{index_inA(j)},"Low Throughput")
%                num_of_low(i)=num_of_low(i)+1;
%            else
%                num_of_both(i)=num_of_both(i)+1;
%  
%            end 
%        end 
%        for j=1:length(index_inB)
%            if strcmp(Interaction_Throughput{index_inB(j)},"High Throughput")
%                num_of_high(i)=num_of_high(i)+1;
%            elseif strcmp(Interaction_Throughput{index_inB(j)},"Low Throughput")
%                num_of_low(i)=num_of_low(i)+1;
%            else
%                num_of_both(i)=num_of_both(i)+1;
%  
%            end 
%        end 
% end
% 
% 
% 
