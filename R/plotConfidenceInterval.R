# The graphical representation of "chipDynoActTransFact's" TF activity! 


plotConfidenceInterval <- function(tf, err){

df <- data.frame(x = 1 : length(tf),
                 F =tf,
                 L =(tf-err),
                 U =(tf+err))

 plot(df$x, df$F, col='red', type = "l")
 #make polygon where coordinates start with lower limit and 
 # then upper limit in reverse order
 polygon(c(df$x,rev(df$x)),c(df$L,rev(df$U)),col = "grey75", border = FALSE)
 lines(df$x, df$F, lwd = 2)
 #add red lines on borders of polygon
 lines(df$x, df$U, col="red",lty=2)
 lines(df$x, df$L, col="red",lty=2)

}

##load("test130313_2.RData")
# source("plotErrorBar.R")
#plotErrorBar(TF[1,],TFError[1,]);

M <- matrix(c(rep(1:24)), byrow=TRUE, nrow=4) # Choose the position by matrix setting!
layout(M) 
for (i in 1:24){
plotConfidenceInterval(TF[i,],TFError[i,])
}

