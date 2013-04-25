x <- rnorm(100) 
M <- matrix(c(rep(1:5, e=2), 6, 7), byrow=TRUE, nrow=2) # Choose the position by matrix setting!
M <- matrix(c(1,2,3,1,4,4), byrow=TRUE, nrow=2) 
M <- matrix(c(1,2,3,1,4,3), byrow=TRUE, nrow=2)

layout(M) 
plot(x) 
hist(x) 
qqnorm(x) 
boxplot(x) 
plot(density(x)) 
plot(abs(x)) 
hist(abs(x)) 
