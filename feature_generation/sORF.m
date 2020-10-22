clc;
clear all;
close all;
a=readtable('full_table.csv');
Seq=table2array(a(1:end,3));
ORF=table2array(a(1:end,2));
Num_of_shifted_ORF=zeros(size(Seq,1),1);
Max_length=zeros(size(Seq,1),1);
Mean_length=zeros(size(Seq,1),1);

for i=1:size(Seq,1)
  seq=Seq{i};
  seq=seq(4:end);
  f=strfind(seq,'ATG');
  Len2=[];
  for j=1:length(f)
     s=seq(f(j):end); 
     Len=0;
     for k=4:3:length(s)-2
         if (s(k)=='T'&&s(k+1)=='G'&&s(k+2)=='A') || (s(k)=='T'&&s(k+1)=='A'&&s(k+2)=='G') || (s(k)=='T'&&s(k+1)=='A'&&s(k+2)=='A')
             stop=k;
             Len=k+2;
             break
         end 
     end
     if Len~=0
         Len2=[Len2,Len];
     end
  end
  if ~isempty(Len2)
      Mean_length(i,1)=mean(Len2);
      Max_length(i,1)=max(Len2);
      Num_of_shifted_ORF(i,1)=length(Len2);
  end
end

Num_of_shifted_ORF_30=zeros(size(Seq,1),1);
Max_length_30=zeros(size(Seq,1),1);
Mean_length_30=zeros(size(Seq,1),1);

%first 30 codons
for i=1:size(Seq,1)
  seq=Seq{i};
  if length(seq)>=90
      seq=seq(4:90);
      f=strfind(seq,'ATG');
      Len2=[];
      for j=1:length(f)
          s=seq(f(j):end); 
          Len=0;
          for k=4:3:length(s)-2
              if (s(k)=='T'&&s(k+1)=='G'&&s(k+2)=='A') || (s(k)=='T'&&s(k+1)=='A'&&s(k+2)=='G') || (s(k)=='T'&&s(k+1)=='A'&&s(k+2)=='A')
                  stop=k;
                  Len=k+2;
                  break
              end
          end
          if Len~=0
              Len2=[Len2,Len];
          end
      end
      if ~isempty(Len2)
          Mean_length_30(i,1)=mean(Len2);
          Max_length_30(i,1)=max(Len2);
          Num_of_shifted_ORF_30(i,1)=length(Len2);
      end
  end
end

Num_of_shifted_ORF_200=zeros(size(Seq,1),1);
Max_length_200=zeros(size(Seq,1),1);
Mean_length_200=zeros(size(Seq,1),1);

%first 200 codons
for i=1:size(Seq,1)
  seq=Seq{i};
  if length(seq)>=600
      seq=seq(4:600);
      f=strfind(seq,'ATG');
      Len2=[];
      for j=1:length(f)
          s=seq(f(j):end); 
          Len=0;
          for k=4:3:length(s)-2
              if (s(k)=='T'&&s(k+1)=='G'&&s(k+2)=='A') || (s(k)=='T'&&s(k+1)=='A'&&s(k+2)=='G') || (s(k)=='T'&&s(k+1)=='A'&&s(k+2)=='A')
                  stop=k;
                  Len=k+2;
                  break
              end
          end
          if Len~=0
              Len2=[Len2,Len];
          end
      end
      if ~isempty(Len2)
          Mean_length_200(i,1)=mean(Len2);
          Max_length_200(i,1)=max(Len2);
          Num_of_shifted_ORF_200(i,1)=length(Len2);
      end
  end
end

final_table1=table(ORF,Mean_length,Max_length,Num_of_shifted_ORF,Mean_length_30,Max_length_30,Num_of_shifted_ORF_30,Mean_length_200,Max_length_200,Num_of_shifted_ORF_200);
writetable(final_table1,'sORF.csv');

%Normalized
Mean_length=(Mean_length-mean(Mean_length))./std(Mean_length);
Max_length=(Max_length-mean(Max_length))./std(Max_length);
Num_of_shifted_ORF=(Num_of_shifted_ORF-mean(Num_of_shifted_ORF))./std(Num_of_shifted_ORF);
Mean_length_30=(Mean_length_30-mean(Mean_length_30))./std(Mean_length_30);
Max_length_30=(Max_length_30-mean(Max_length_30))./std(Max_length_30);
Num_of_shifted_ORF_30=(Num_of_shifted_ORF_30-mean(Num_of_shifted_ORF_30))./std(Num_of_shifted_ORF_30);
Mean_length_200=(Mean_length_200-mean(Mean_length_200))./std(Mean_length_200);
Max_length_200=(Max_length_200-mean(Max_length_200))./std(Max_length_200);
Num_of_shifted_ORF_200=(Num_of_shifted_ORF_200-mean(Num_of_shifted_ORF_200))./std(Num_of_shifted_ORF_200);

final_table2=table(ORF,Mean_length,Max_length,Num_of_shifted_ORF,Mean_length_30,Max_length_30,Num_of_shifted_ORF_30,Mean_length_200,Max_length_200,Num_of_shifted_ORF_200);
writetable(final_table2,'sORF_Normalized.csv');
