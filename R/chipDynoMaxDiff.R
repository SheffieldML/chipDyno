# CHIPDYNOMAXDIFF computes most significant changes in TFAs
# CHIPDYNO toolbox
# chipDynoMaxDiff.R version 1.0.1
# FORMAT chipDynoMaxDiff <- function(TF,TFErrorDiff)
# DESC computes most significant changes in TFAs
# ARG TF: gene specific transcription factor activity
# ARG TFErrorDiff: error in gene specific transcription factor activity
# RETURN f : most significant changes in TFAs
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipDynoExpectationsFast, chipDynoExpectationsFastNoise

chipDynoMaxDiff <- function(TF,TFErrorDiff) {

nTargets=nrow(TF);
npts=ncol(TF);
preDiffs=array(0, dim <- c(npts,npts));
diffs=array(0, dim <- c(npts,npts,nTargets));
f= array(0, dim <- c(1,nTargets));

for (i in 1: nTargets) {
	for (j in 2:(npts-1)) {
		for (l in j : npts) {
			preDiffs[j,l]=TF[i,j]-TF[i,l];
		}
	}
	diffs[ , ,i] = preDiffs - t(preDiffs);
	f[i]=max(diffs[,,i]/TFErrorDiff[,,i]);  ## ?? ./
}

return(f)
}