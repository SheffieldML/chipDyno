# Load the library
library("foreach",lib.loc="/home/muhammad/R/installed_packages/")
library("iterators",lib.loc="/home/muhammad/R/installed_packages/")
library("snow",lib.loc="/home/muhammad/R/installed_packages/")
library("doSNOW",lib.loc="/home/muhammad/R/installed_packages/")

#library(doSNOW)
# rm(list=ls())

cl<-makeCluster(4) # I have two cores
registerDoSNOW(cl)

table <- data.frame(a=rnorm(10),b=rnorm(10))

process <- function(table)
              {for (loop in (1:nrow(table)))
                   {#table[loop,"c"] <- with(table[loop,], a*b)
			table[loop,"c"] <- (table[loop,1]+table[loop,2])
                    assign("table",table,envir=.GlobalEnv)
                   }
              }

system.time(process(table))

system.time(foreach(j=1:2 ) %dopar% process(table))

stopCluster(cl)
