clc
clear all
mergedFiles = dir('*.merged'); 
numfiles = length(mergedFiles);
mydata = cell(1, numfiles);
columns_with_scores=[3,4,6,12];
column_with_prediction=[7,8];
gene_names=cell(1, numfiles);
mean_score=zeros(5839,1);%mean on all the columns that has numeric values- 3,4,6,12
disorder_percentage1=zeros(5839,1);%first 4 algorithm
disorder_percentage2=zeros(5839,1);%next 3 algorithms
disorder_percentage_both_algorithms=zeros(5839,1);%percentage 1 and 2 were disordered
anchor_percentage=zeros(5839,1);%column 13
anchor_disorder_percentage=zeros(5839,1);%column 12
IUPRED_percentage=zeros(5839,1);%column 3
VSL2B_percentage=zeros(5839,1);%column 4
DisEMBL_percentage=zeros(5839,1);%column 5
MoreRONN_percentage=zeros(5839,1);%column 6
disorder30_percentage=zeros(5839,1);%column 9

for k =1:numfiles 
    fileid=fopen(mergedFiles(k).name);
    file_name=split(mergedFiles(k).name,"_"); 
    gene_names{k}=file_name{1};
    mydata{k}=textscan(fileid, '%f %s %f %f %s %f %f %s %s %f %s %f %f','HeaderLines',1,'Delimiter','\t','EndOfLine','\n');
    mean_score(k)=mean([mean(mydata{k}{1,3}),mean(mydata{k}{1,4}),mean(mydata{1,k}{1,6}),mean(mydata{1,k}{1,12})]);
    IUPRED_percentage(k)=mean(mydata{1,k}{1,3});
    VSL2B_percentage(k)=mean(mydata{1,k}{1,4});
    MoreRONN_percentage(k)=mean(mydata{1,k}{1,6});
    anchor_percentage(k)=mean(mydata{1,k}{1,12});
    for i=1:length( mydata{k}{1})
        if (mydata{k}{7}(i)==3 | mydata{k}{7}(i)==4) && string(mydata{k}{8}(i))=='DISORDERED'
            disorder_percentage_both_algorithms(k)=disorder_percentage_both_algorithms(k)+1;
        end 
        if mydata{k}{7}(i)==3 | mydata{k}{7}(i)==4
            disorder_percentage1(k)=disorder_percentage1(k)+1;
        end 
        if string(mydata{k}{8}(i))=='DISORDERED'
            disorder_percentage2(k)=disorder_percentage2(k)+1;
        end
        if mydata{k}{13}(i)==1
            anchor_disorder_percentage(k)=anchor_disorder_percentage(k)+1;
        end
        if string(mydata{k}{5}(i))=='DISORDERED'
            DisEMBL_percentage(k)=DisEMBL_percentage(k)+1;
        end
        if string(mydata{k}{9}(i))=='DIS_30'
            disorder30_percentage(k)=disorder30_percentage(k)+1;
        end
    end

        disorder_percentage_both_algorithms(k)=disorder_percentage_both_algorithms(k)/length(mydata{k}{1});
        disorder_percentage1(k)=disorder_percentage1(k)/length(mydata{k}{1});
        disorder_percentage2(k)=disorder_percentage2(k)/length(mydata{k}{1});
        anchor_disorder_percentage(k)=anchor_disorder_percentage(k)/length(mydata{k}{1});
        DisEMBL_percentage(k)=DisEMBL_percentage(k)/length(mydata{k}{1});
        disorder30_percentage(k)=disorder30_percentage(k)/length(mydata{k}{1});
        fclose(fileid);
end

%first 30 windows, each window with length 50
disorder_percentage1_first30=zeros(5839,30);
disorder_std1_first30=zeros(5839,30);
disorder_percentage2_first30=zeros(5839,30);
disorder_std2_first30=zeros(5839,30);
disorder30_percentage_first30=zeros(5839,30);
disorder30_std_first30=zeros(5839,30);
disorder_anchor_percentage_first30=zeros(5839,30);
disorder_anchor_std_first30=zeros(5839,30);
max_disorder_percentage_first30=zeros(5839,30);

for k=1:numfiles
    if length(mydata{k}{1})<79
        max_i=length(mydata{k}{1})-50;
        disorder_percentage1_first30(k,max_i+1:30)=NaN;
        disorder_std1_first30(k,max_i+1:30)=NaN;
        disorder_percentage2_first30(k,max_i+1:30)=NaN;
        disorder_std2_first30(k,max_i+1:30)=NaN;
        disorder30_percentage_first30(k,max_i+1:30)=NaN;
        disorder30_std_first30(k,max_i+1:30)=NaN;
        disorder_anchor_percentage_first30(k,max_i+1:30)=NaN;
        disorder_anchor_std_first30(k,max_i+1:30)=NaN;
        max_disorder_percentage_first30(k,max_i+1:30)=NaN;
    else
        max_i=30;
    end
    for i=1:max_i
        disorder_percentage1_first30(k,i)=sum(mydata{k}{7}(i:i+49)==3)+sum(mydata{k}{7}(i:i+49)==4);
        disorder_std1_first30(k,i)=std([mydata{k}{7}(i:i+49)==3]+[mydata{k}{7}(i:i+49)==4]);
        disorder_percentage2_first30(k,i)=sum(string(mydata{k}{8}(i:i+49))=="DISORDERED");
        disorder_std2_first30(k,i)=std(string(mydata{k}{8}(i:i+49))=="DISORDERED");
        disorder30_percentage_first30(k,i)=sum(string(mydata{k}{9}(i:i+49))=="DIS_30");
        disorder30_std_first30(k,i)=std(string(mydata{k}{9}(i:i+49))=="DIS_30");
        disorder_anchor_percentage_first30(k,i)=sum(mydata{k}{13}(i:i+49));
        disorder_anchor_std_first30(k,i)=std(mydata{k}{13}(i:i+49));
        max_disorder_percentage_first30(k,i)=max([max(mydata{1,k}{1,3}(i:i+49)) max(mydata{1,k}{1,4}(i:i+49)) max(mydata{1,k}{1,6}(i:i+49)) max(mydata{1,k}{1,12}(i:i+49))]);   
    end
      disorder_percentage1_first30(k,1:30)=disorder_percentage1_first30(k,1:30)/50;
      disorder_percentage2_first30(k,1:30)=disorder_percentage2_first30(k,1:30)/50;
      disorder30_percentage_first30(k,1:30)=disorder30_percentage_first30(k,1:30)/50;
      disorder_anchor_percentage_first30(k,1:30)=disorder_anchor_percentage_first30(k,1:30)/50;
end

disorder_percentage1_lastwindow=zeros(5839,1);
disorder_std1_lastwindow=zeros(5839,1);
disorder_percentage2_lastwindow=zeros(5839,1);
disorder_std2_lastwindow=zeros(5839,1);
disorder30_percentage_lastwindow=zeros(5839,1);
disorder30_std_lastwindow=zeros(5839,1);
disorder_anchor_percentage_lastwindow=zeros(5839,1);
disorder_anchor_std_lastwindow=zeros(5839,1);
max_disorder_percentage_lastwindow=zeros(5839,1);

for k=1:numfiles
    length_of_file=length(mydata{k}{1});
    disorder_percentage1_lastwindow(k)=(sum(mydata{k}{7}(length_of_file-49:length_of_file)==3)+sum(mydata{k}{7}(length_of_file-49:length_of_file)==4))/50;
    disorder_percentage1_lastwindow(k)=std([mydata{k}{7}(length_of_file-49:length_of_file)==3]+[mydata{k}{7}(length_of_file-49:length_of_file)==4]);
    disorder_percentage2_lastwindow(k)=(sum(string(mydata{k}{8}(length_of_file-49:length_of_file))=="DISORDERED"))/50;
    disorder_std2_lastwindow(k)=std(string(mydata{k}{8}(length_of_file-49:length_of_file))=="DISORDERED");
    disorder30_percentage_lastwindow(k)=(sum(string(mydata{k}{9}(length_of_file-49:length_of_file))=="DIS_30"))/50;
    disorder30_std_lastwindow(k)=std(string(mydata{k}{9}(length_of_file-49:length_of_file))=="DIS_30");
    disorder_anchor_percentage_lastwindow(k)=(sum(mydata{k}{13}(length_of_file-49:length_of_file)))/50;
    disorder_anchor_std_lastwindow(k)=std(mydata{k}{13}(length_of_file-49:length_of_file));
    max_disorder_percentage_lastwindow(k)=max([max(mydata{1,k}{1,3}(length_of_file-49:length_of_file)),max(mydata{1,k}{1,4}(length_of_file-49:length_of_file)),max(mydata{1,k}{1,6}(length_of_file-49:length_of_file)),max(mydata{1,k}{1,12}(length_of_file-49:length_of_file))]);   
end

        