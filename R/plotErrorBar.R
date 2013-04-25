# The graphical representation of "chipDynoActTransFact's" TF activity! 

plotErrorBar <- function(y,SE){

add.error.bars <- function(x,y,SE,w,col=1){
x0 = x; y0 = (y-SE); x1 =x; y1 = (y+SE);
arrows(x0, y0, x1, y1, code=3,angle= 90,length=w,col=col);
}

### Test
x <- c(1:length(y));
#Y <- TF[1,];
#SE <- TFError[1,];
plot(x,y, type = 'l', col='green4', las=1);
add.error.bars(x,y,SE,0.05,col='red');

###
}

##load("test130313_2.RData")
# source("plotErrorBar.R")
#plotErrorBar(TF[1,],TFError[1,]);

M <- matrix(c(rep(1:24)), byrow=TRUE, nrow=4) # Choose the position by matrix setting!
layout(M) 
for (i in 1:24){
plotErrorBar(TF[i,],TFError[i,])
}

