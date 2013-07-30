#function f=chipDynoMaxDiff(TF,TFErrorDiff);

#% CHIPDYNOMAXDIFF computes most significant changes in TFAs
#%
#%	Description:
#%	f=chipDynoMaxDiff(TF,TFErrorDiff);
#%% 	chipDynoMaxDiff.R version 0.1.0
##

chipDynoMaxDiff = function(TF,TFErrorDiff) {

nTargets=nrow(TF);
npts=ncol(TF);
preDiffs=array(0, dim <- c(npts,npts));
diffs=array(0, dim <- c(npts,npts,nTargets));
f= array(0, dim <- c(1,nTargets));

#nTargets=size(TF,1);
#npts=size(TF,2);
#preDiffs=zeros(npts,npts);
#diffs=zeros(npts,npts,nTargets);
#f=zeros(1,nTargets);


for (i in 1: nTargets) {
	for (j in 2:(npts-1)) {
		for (l in j : npts) {
			preDiffs[j,l]=TF[i,j]-TF[i,l];
		}
	}
	diffs[ , ,i] = preDiffs - t(preDiffs);
	f[i]=max(diffs[,,i]/TFErrorDiff[,,i]);  ## ?? ./
}



#for i=1:nTargets
#    for j=2:npts-1
#        for l=j:npts
#            preDiffs(j,l)=TF(i,j)-TF(i,l);
#        end
#    end
#    diffs(:,:,i)=preDiffs-preDiffs';
#    f(i)=max(max(diffs(:,:,i)./TFErrorDiff(:,:,i)));
#end

return(f)
}
