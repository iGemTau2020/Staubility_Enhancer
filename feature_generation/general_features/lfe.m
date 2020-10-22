clc;
clear all;

%ORF=csvread('export_taxid_559292_begin_dlfe.csv',1,0,[1,0,5797,0]);
info_dlfe=csvread('export_taxid_559292_begin_dlfe.csv',1,1,[1,1,5797,31]);
info_native=csvread('export_taxid_559292_begin_native.csv',1,1,[1,1,5797,31]);
info_shuffled=csvread('export_taxid_559292_begin_shuffled.csv',1,1,[1,1,5797,31]);
info_windows=csvread('export_taxid_559292_begin_shuffled.csv',0,1,[0,1,0,31]);
info_dlfe=info_dlfe./max(abs(info_dlfe));
info_native=info_native./max(abs(info_native));
info_shuffled=info_shuffled./max(abs(info_shuffled));
mean_dlfe=mean(info_dlfe,2);
mean_native=mean(info_native,2);
mean_shuffled=mean(info_shuffled,2);
a=readtable('export_taxid_559292_begin_dlfe.csv');
ORF=a(1:end,1);
ORF=table2cell(ORF);
ORF{1}='ORF';
ORF2=table2cell(a(2:end,1));
final_table1=table(ORF,[info_windows;info_dlfe]);
final_table2=table(ORF,[info_windows;info_native]);
final_table3=table(ORF,[info_windows;info_shuffled]);
final_table4=table(ORF2,mean_native,mean_shuffled,mean_dlfe,'VariableNames',{'ORF','mean_native','mean_shuffled','mean_dlfe'});
writetable(final_table1,'dlfe.csv');
writetable(final_table2,'lfe_native.csv');
writetable(final_table3,'lfe_shuffled.csv');
writetable(final_table4,'mean_lfe.csv');

