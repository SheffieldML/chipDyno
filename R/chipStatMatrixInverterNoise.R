# CHIPSTATMATRIXINVERTERNOISE inverts block tridiagonal matrices for chipChip
# CHIPDYNO toolbox
# chipStatMatrixInverterNoise.R version 1.0.1
# FORMAT chipStatMatrixInverterNoise <- function(Sigma, gamma, beta, precs, x, npts)
# DESC inverts block tridiagonal matrices for chipChip
# ARG Sigma : prior covariance matrix of TFA
# ARG gamma : degree of temporal continuity
# ARG beta :
# ARG precs : uncertainty of the expression level
# ARG x : connectivity measurement between genes and transcription factors
# ARG npts : number of transcription factors
# RETURN f : inverted block tridiagonal matrices 
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipStatMatrixInverter

chipStatMatrixInverterNoise <- function(Sigma, gamma, beta, precs, x, npts){

lambda = t(x) %*% Sigma %*% x
lambda = as.vector(lambda)
nTrans = length(x)
Y = Sigma %*% x
factor=cos(gamma)^2
UcoeffInvSigma=mat.or.vec(1,npts)
UcoeffXXT=mat.or.vec(1,npts)
LcoeffId=mat.or.vec(1,npts-1) 		# computes LU dec exploiting simple block
LcoeffXYT=mat.or.vec(1,npts-1)

UcoeffXXT[1]=(beta^-2+precs[1]^-1)^-1
UcoeffInvSigma[1]=(1-factor^2)^-1
LcoeffXYT[1]=factor*(1-factor^2)^-1*(beta^-2+precs[1]^-1)^-1/
	(UcoeffInvSigma[1]* (UcoeffInvSigma[1]+(beta^-2+precs[1]^-1)^-1*lambda));
LcoeffId[1]=-factor

for (i in 2:(npts-1)){
	UcoeffXXT[i]=(beta^-2+precs[i]^-1)^-1+factor*(1-factor^2)^-1*LcoeffXYT[(i-1)]
	UcoeffInvSigma[i]=(1+factor^2)*(1-factor^2)^-1+
			factor*(1-factor^2)^-1*LcoeffId[(i-1)]
	LcoeffXYT[i]=factor*(1-factor^2)^-1*UcoeffXXT[i]/
			(UcoeffInvSigma[i]*(UcoeffInvSigma[i]+UcoeffXXT[i]*lambda))
	LcoeffId[i]=-factor*(1-factor^2)^-1*UcoeffInvSigma[i]^-1
}

UcoeffXXT[ncol(UcoeffXXT)]=(beta^-2+precs[length(precs)]^-1)^-1+
			factor*(1-factor^2)^-1*LcoeffXYT[ncol(LcoeffXYT)]
UcoeffInvSigma[ncol(UcoeffInvSigma)]=(1-factor^2)^-1+
			factor*(1-factor^2)^-1*LcoeffId[ncol(LcoeffId)]
#%lambda=Y'*x;
invL.XYT=mat.or.vec(npts,npts); #%computes the inverse of the L bit
invL.Id=mat.or.vec(npts,npts)

for (i in 1:npts){
	invL.Id[i,i]=1
}

for (i in 2:npts) {
	invL.Id[i,(i-1)]=(-1)^(2*i-1)*LcoeffId[i-1]
	invL.XYT[i,(i-1)]=(-1)^(2*i-1)*LcoeffXYT[i-1]
}

for (i in 3 : npts) {
	for (j in 1:(i-2)){
		invL.Id[i,j]=invL.Id[(i-1),j]*invL.Id[i,(i-1)]
		invL.XYT[i,j]=(invL.Id[(i-1),j]*invL.XYT[i,(i-1)]+ 
			invL.Id[i,(i-1)]*invL.XYT[(i-1),j]+ 
			invL.XYT[i,(i-1)]*invL.XYT[(i-1),j]*lambda)
	}
}

invU.Sigma=mat.or.vec(npts,npts)
invU.YYT=mat.or.vec(npts,npts)
for (i in 1:(npts-1)){
	invU.Sigma[i,i]=-LcoeffId[i]*(1-factor^2)/factor
	invU.YYT[i,i]=-LcoeffXYT[i]*(1-factor^2)/factor
}

invU.Sigma[nrow(invU.Sigma),ncol(invU.Sigma)]=
	UcoeffInvSigma[length(UcoeffInvSigma)]^-1

invU.YYT[nrow(invU.YYT),ncol(invU.YYT)]=-UcoeffXXT[length(UcoeffXXT)]/
	(UcoeffInvSigma[length(UcoeffInvSigma)]*
			(UcoeffInvSigma[length(UcoeffInvSigma)]+
				UcoeffXXT[length(UcoeffXXT)]*lambda))

for (i in 1:(npts-1)){
	for (j in 1:i){
		invU.Sigma[(npts-i),(npts-j+1)]=factor*(1-factor^2)^-1* 
			invU.Sigma[(npts-i+1),(npts-j+1)]*invU.Sigma[(npts-i),(npts-i)];
        invU.YYT[(npts-i),(npts-j+1)]=factor*(1-factor^2)^-1*
        	(invU.Sigma[(npts-i+1),(npts-j+1)]*invU.YYT[(npts-i),(npts-i)]+
        		invU.Sigma[(npts-i),(npts-i)]*invU.YYT[(npts-i+1),(npts-j+1)]+ 
        			invU.YYT[(npts-i+1),(npts-j+1)]*
        				invU.YYT[(npts-i),(npts-i)]*lambda)
	}
}

invC.Sigma=mat.or.vec(npts,npts); #%computes the inverses of the matrix;
invC.YYT=mat.or.vec(npts,npts)

for (i in 1:npts) {
	for (j in 1:npts) {
        	invC.Sigma[i,j]=invU.Sigma[i,]%*%invL.Id[,j]
		invC.YYT[i,j]=invU.Sigma[i,]%*%invL.XYT[,j]+invU.YYT[i,]%*% invL.Id[,j]+
			lambda*invU.YYT[i,]%*%invL.XYT[,j]
	}
}

invC <- list(invC.Sigma, invC.YYT)

return(invC)
}