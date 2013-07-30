#Sample#1
rm(list=ls())

require("Matrix")

tfId=8
#tfId=33 # To examin the activity of TransNames[tfId,]

source("chipDynoTransFactNoise.R")
source("chipDynoMaxDiff.R")

load("/home/muhammad/H-drive/CElegans/Results_cElegans_Optim_Sample1.RData")

nTrans=nrow(TransNames);
lst=list();
newX=array(0, dim <-c(dim(X)));
newXVals=array(0, dim <-c(dim(X)));

i=tfId
expectations =chipDynoTransFactNoise(data,X,Sigma,beta, precs, gamma,mu, TransNames, annotation, TransNames[i,]);
S1_TF = expectations[[1]]
S1_TFError = expectations [[2]]
S1_TFErrorDiff = expectations [[3]]

###################

#Sample#2
load("/home/muhammad/H-drive/CElegans/Results_cElegans_Optim_Sample2.RData")

nTrans=nrow(TransNames);
lst=list();
newX=array(0, dim <-c(dim(X)));
newXVals=array(0, dim <-c(dim(X)));

i=tfId
expectations =chipDynoTransFactNoise(data,X,Sigma,beta, precs, gamma,mu, TransNames, annotation, TransNames[i,]);
S2_TF = expectations[[1]]
S2_TFError = expectations [[2]]
S2_TFErrorDiff = expectations [[3]]

############

#Sample#3
load("/home/muhammad/H-drive/CElegans/Results_cElegans_Optim_Sample3.RData")

nTrans=nrow(TransNames);
lst=list();
newX=array(0, dim <-c(dim(X)));
newXVals=array(0, dim <-c(dim(X)));

i=tfId
expectations =chipDynoTransFactNoise(data,X,Sigma,beta, precs, gamma,mu, TransNames, annotation, TransNames[i,]);
S3_TF = expectations[[1]]
S3_TFError = expectations [[2]]
S3_TFErrorDiff = expectations [[3]]

###### Find the geneNames####

i=tfId
name = TransNames[i,]
index=which(name == TransNames);
genesIn=which(X[,index]!=0);
actGeneNames = annotation[genesIn]

####

# The graphical representation of "chipDynoActTransFact's" TF activity! 

plotErrorBar <- function(y,SE,z,SE_z,zz, SE_zz, actGeneName){

add.error.bars <- function(x,y,SE,w,col=1){
x0 = x; y0 = (y-SE); x1 =x; y1 = (y+SE);
arrows(x0, y0, x1, y1, code=3,angle= 90,length=w,col=col);
}

### Test
x <- c(1:length(y));
#x <- c(3,20,24,48,72);
#Y <- TF[1,];
#SE <- TFError[1,];
#plot(x,y, type = 'l', col='green4', las=1);
#plot(x,y, type = 'l', col='green4', las=1, main/sub="gene specific TFA", xlab="t", ylab="TFA");

plot(x,y, type = 'l', col='green4', las=1, ylim=c(0,13.5), xlab="t", ylab="TFA", main=bquote("Gene : " ~ .(actGeneName)));
#plot(x,y, type = 'l', col='green4', las=1, ylim=c(0,13.5), xlab="t", ylab="TFA", sub=bquote("Gene Specific TFA of T20B12.8"), main=bquote("Gene Name = " ~ .(actGeneName)));
#plot(x,y, type = 'l', col='green4', las=1, ylim=c(0,13.5), xlab="t", ylab="TFA"));

lines(x,z, type = 'l', col='blue', las=1);
lines(x,zz, type = 'l', col='black', las=1);

#lnames <- c('Exp1', 'Exp2', 'Exp3') 
legend('topright', legend <- c('Exp1', 'Exp2', 'Exp3'), col=c("green4","blue", "black"), lty = c('solid', 'solid', 'solid')) 

add.error.bars(x,y,SE,0.05,col='red');
add.error.bars(x,z,SE_z,0.05,col='red');
add.error.bars(x,zz,SE_zz,0.05,col='red');

}

##load("test130313_2.RData")
# source("plotErrorBar.R")
#plotErrorBar(TF[1,],TFError[1,]);

M <- matrix(c(rep(1:6)), byrow=TRUE, nrow=2) # Choose the position by matrix setting!
layout(M) 

#jpeg('rplot.jpg')
#bmp("testplot.bmp")
#png(filename="TFA_data_sample1.png")
for (i in 1:6){
plotErrorBar(S1_TF[i,],S1_TFError[i,],S2_TF[i,],S2_TFError[i,],S3_TF[i,],S3_TFError[i,], actGeneNames[i])
}
#dev.off()

# setEPS()
# postscript("T20B12_8_1.eps")
# source("plotErrorBar_Multiple.R")
# dev.off()