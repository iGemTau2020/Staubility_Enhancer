function CFS=calculate_CFS(rcf,rff,feature_ind)
% rcf = feature-class correlation matrix
% rff = feature-feature correlation matrix
partial_rcf=abs(rcf(feature_ind));
partial_rff=abs(rff(feature_ind,feature_ind));
CFS=sum(partial_rcf)/sqrt(sum(sum(partial_rff)));
