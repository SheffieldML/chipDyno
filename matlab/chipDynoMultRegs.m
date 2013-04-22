function list=chipDynoMultRegs(lista,data,X,Sigma,beta,precs,gamma, ...
                              mu,TransNames,annotation);

% CHIPDYNOMULTREGS finds genes that have multiple regulators
%
%	Description:
%	list=chipDynoMultRegs(lista,data,X,Sigma,beta,precs,gamma, ...
%                              mu,TransNames,annotation);
%% 	chipDynoMultRegs.m version 1.4

list=[];
for i=1:size(lista,1)
  [TF,TFError]=chipDynoTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                      TransNames, annotation, lista(i));
  index=find(strcmp(lista(i),TransNames));
  anno=annotation(find(X(:,index)));
  [values,indices]=sort(sqrt(var(TF'))./mean(TFError,2)','descend');
  anns=anno(indices(find(values>2)));
  lister=chipDynoPerGenes(anns);
   counter=0;
  for j=1:size(lister,1)
   
    if strcmp(lister(j),'')==0
      counter=counter+1;
       if sum(X(find(strcmp(lister(j),annotation)),:))>1
         list=[list;lister(j)];
       end
    end
    
  end
end

    
    
    
    
