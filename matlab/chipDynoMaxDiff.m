function f=chipDynoMaxDiff(TF,TFErrorDiff);

% CHIPDYNOMAXDIFF computes most significant changes in TFAs
%
%	Description:
%	f=chipDynoMaxDiff(TF,TFErrorDiff);
%% 	chipDynoMaxDiff.m version 1.2


nTargets=size(TF,1);
npts=size(TF,2);
preDiffs=zeros(npts,npts);
diffs=zeros(npts,npts,nTargets);
f=zeros(1,nTargets);
for i=1:nTargets
    for j=2:npts-1
        for l=j:npts
            preDiffs(j,l)=TF(i,j)-TF(i,l);
        end
    end
    diffs(:,:,i)=preDiffs-preDiffs';
    f(i)=max(max(diffs(:,:,i)./TFErrorDiff(:,:,i)));
end
