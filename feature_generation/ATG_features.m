clc;
clear all;
close all;
a=readtable('full_table.csv');
ORF=table2array(a(1:end,3));
Gene_Name=table2array(a(1:end,2));

%Number of lt ATG ORF
Numer_of_alt_ATG_ORF=zeros(length(Gene_Name),1);
Numer_of_alt_ATG_ORF_30codons=zeros(length(Gene_Name),1);
Numer_of_alt_ATG_ORF_200codons=zeros(length(Gene_Name),1);

for i=1:size(ORF,1)
  O=ORF{i};
  O=O(4:end);
  Numer_of_alt_ATG_ORF(i,1)=length(strfind(O,'ATG'));
  if length(O)>=600
      Numer_of_alt_ATG_ORF_200codons(i,1)=length(strfind(O(1:200),'ATG'));
  end
  if length(O)>=90
      Numer_of_alt_ATG_ORF_30codons(i,1)=length(strfind(O(1:90),'ATG'));
  end
end


Numer_of_alt_ATG_ORF_30codons(isnan(Numer_of_alt_ATG_ORF_30codons))=0;
Numer_of_alt_ATG_ORF_200codons(isnan(Numer_of_alt_ATG_ORF_200codons))=0;
Numer_of_alt_ATG_ORF(isnan(Numer_of_alt_ATG_ORF))=0;

a=readtable('ATG_Tamir.xlsx');
gene_name_ORF=table2array(a(7139:end,2));
gene_name_5UTR=table2array(a(1:1277,2));
Gene_Name=table2array(a(1278:7138,2));
Gene_Location_5UTR=table2array(a(1:1277,3));
Gene_Location_ORF=table2array(a(7139:end,3));

Absolute_Context_Score_Alt_5UTR=exp(table2array(a(1:1277,4)));
Absolute_Context_Score_Alt_ORF=exp(table2array(a(7139:end,4)));
Relative_Context_Score_Alt_5UTR=exp(table2array(a(1:1277,5)));
Relative_Context_Score_Alt_ORF=exp(table2array(a(7139:end,5)));

%ATG score for Main ATG
Absolute_Context_Score_Main_ATG=exp(table2array(a(1278:7138,4)));

%Max/Mean absolute score
Max_ATG_Context_Score_ORF_30codons=zeros(length(Gene_Name),1);
Mean_ATG_Context_Score_ORF_30codons=zeros(length(Gene_Name),1);

%Max/Mean relative score
Max_relative_ATG_Context_Score_ORF_30codons=zeros(length(Gene_Name),1);
Mean_relative_ATG_Context_Score_ORF_30codons=zeros(length(Gene_Name),1);

%Number of alt ATG 5UTR
Numer_of_alt_ATG_5UTR=zeros(length(Gene_Name),1);
Numer_of_alt_ATG_5UTR_30codons=zeros(length(Gene_Name),1);

%Distance of 1st ATG in the transcript segment from the main START ATG in the ORF
Dist_5UTR_to_Main_ATG=zeros(length(Gene_Name),1);

for i=1:size(Gene_Name,1)
    Numer_of_alt_ATG_5UTR(i,1)=length(find(contains(gene_name_5UTR,Gene_Name(i))));
    loc_ORF=Gene_Location_ORF(find(contains(gene_name_ORF,Gene_Name(i))));
    val_ORF_a=Absolute_Context_Score_Alt_ORF(find(contains(gene_name_ORF,Gene_Name(i))));
    val_ORF_r=Relative_Context_Score_Alt_ORF(find(contains(gene_name_ORF,Gene_Name(i))));
    if ~isempty(val_ORF_a)
        if ~isempty(find(loc_ORF<=90))
            Max_ATG_Context_Score_ORF_30codons(i,1)=max(val_ORF_a(find(loc_ORF<=90)));
            Mean_ATG_Context_Score_ORF_30codons(i,1)=mean(val_ORF_a(find(loc_ORF<=90)));
            Max_relative_ATG_Context_Score_ORF_30codons(i,1)=max(val_ORF_r(find(loc_ORF<=90)));
            Mean_relative_ATG_Context_Score_ORF_30codons(i,1)=mean(val_ORF_r(find(loc_ORF<=90)));
        end
    end
end


Numer_of_alt_ATG_5UTR(isnan(Numer_of_alt_ATG_5UTR))=0;
MeanMax=[Max_ATG_Context_Score_ORF_30codons,Mean_ATG_Context_Score_ORF_30codons,Max_relative_ATG_Context_Score_ORF_30codons,Mean_relative_ATG_Context_Score_ORF_30codons];

for i=1:size(MeanMax,2)
    b=MeanMax(1:end,i);
    m=min(b(find(b~=0)));
    b(find(b==0))=m/(10^10);
    b=log(b);
    MeanMax(1:end,i)=b;
    mm(i)=log(m);
end

c=readtable('full_table.csv');
ORF=table2array(c(1:end,2));

%fill for all genes
Absolute_Context_Score_Main_ATG=log(Absolute_Context_Score_Main_ATG);
ATG_not_final=[Absolute_Context_Score_Main_ATG,Numer_of_alt_ATG_5UTR,MeanMax];
ATG_final=zeros(size(ORF,1),size(ATG_not_final,2));

%mm(5)=min(Absolute_Context_Score_Main_ATG);
%fill with nan

for i=1:size(ORF,1)
    z=ORF{i};
    check=0;
    if sum(contains(Gene_Name,z))>0
        for j=1:length(find(contains(Gene_Name,z)))
            o=find(contains(Gene_Name,z));
            if strcmp(Gene_Name{o(j)},z)==1
                ATG_final(i,1:end)=ATG_not_final(o(j),1:end);
                check=1;
            end
        end
    else
      ATG_final(i,1:end)=nan;   
    end
    if check==0
        ATG_final(i,1:end)=nan; 
    end
end


Absolute_Context_Score_Main_ATG=ATG_final(1:end,1);
%Absolute_Context_Score_Main_ATG(isnan(Absolute_Context_Score_Main_ATG))=mm(5)/(10^10);
%Absolute_Context_Score_Main_ATG=log(Absolute_Context_Score_Main_ATG);
Numer_of_alt_ATG_5UTR=ATG_final(1:end,2);
Max_ATG_Context_Score_ORF_30codons=ATG_final(1:end,3);
Max_ATG_Context_Score_ORF_30codons(isnan(Max_ATG_Context_Score_ORF_30codons))=mm(1);
Mean_ATG_Context_Score_ORF_30codons=ATG_final(1:end,4);
Mean_ATG_Context_Score_ORF_30codons(isnan(Mean_ATG_Context_Score_ORF_30codons))=mm(2);
Max_relative_ATG_Context_Score_ORF_30codons=ATG_final(1:end,5);
Max_relative_ATG_Context_Score_ORF_30codons(isnan(Max_relative_ATG_Context_Score_ORF_30codons))=mm(3);
Mean_relative_ATG_Context_Score_ORF_30codons=ATG_final(1:end,6);
Mean_relative_ATG_Context_Score_ORF_30codons(isnan(Mean_relative_ATG_Context_Score_ORF_30codons))=mm(4);

Number_of_alt_ATG_tot=Numer_of_alt_ATG_5UTR+Numer_of_alt_ATG_ORF;

final_table1=table(Absolute_Context_Score_Main_ATG,Numer_of_alt_ATG_5UTR,Numer_of_alt_ATG_ORF,Numer_of_alt_ATG_ORF_30codons,Numer_of_alt_ATG_ORF_200codons,Number_of_alt_ATG_tot,Max_ATG_Context_Score_ORF_30codons,Mean_ATG_Context_Score_ORF_30codons,Max_relative_ATG_Context_Score_ORF_30codons,Mean_relative_ATG_Context_Score_ORF_30codons);
writetable(final_table1,'ATG_withnan1.csv');

Numer_of_alt_ATG_5UTR(isnan(Numer_of_alt_ATG_5UTR))=0;
Number_of_alt_ATG_tot=Numer_of_alt_ATG_5UTR+Numer_of_alt_ATG_ORF;

final_table1=table(Absolute_Context_Score_Main_ATG,Numer_of_alt_ATG_5UTR,Numer_of_alt_ATG_ORF,Numer_of_alt_ATG_ORF_30codons,Numer_of_alt_ATG_ORF_200codons,Number_of_alt_ATG_tot,Max_ATG_Context_Score_ORF_30codons,Mean_ATG_Context_Score_ORF_30codons,Max_relative_ATG_Context_Score_ORF_30codons,Mean_relative_ATG_Context_Score_ORF_30codons);
writetable(final_table1,'ATG_withnan.csv');