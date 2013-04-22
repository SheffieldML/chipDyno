function values=chipDynoSignRels(indicator, newX, data, X, Sigma, beta, gamma, mu, TransNames, annotation)

% CHIPDYNOSIGNRELS finds the regulatory strength conserved rels
%
%	Description:
%	values=chipDynoSignRels(indicator, newX, data, X, Sigma, beta, gamma, mu, TransNames, annotation)
%% 	chipDynoSignRels.m version 1.2

[row, col]= find(newX);
values=[];
effectRow=row(find(indicator));
effectCol=col(find(indicator));
for i=1:sum(indicator)
    [TF,TFError]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
                                         TransNames, annotation, ...
                                         TransNames(effectCol(i)),annotation(effectRow(i)));
    values=[values, max(abs((TF-mu(effectCol(i))*ones(size(TF)))./TFError))];
end
