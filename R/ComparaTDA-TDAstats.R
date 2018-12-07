library(TDA)
library(TDAstats)
library(tidyverse)

n <- 200
data <- cbind(runif(n), runif(n), runif(n))

# comparison TDA vs TDAstats
dimension <- 3
xmax <- 1

# TDA
# startTDA <- Sys.time()
# ripsTDA <- ripsDiag(data, 
#                     maxdimension = dimension, 
#                     maxscale = xmax)
# endTDA <- Sys.time()
# print(c('TDA:      ', endTDA-startTDA))
# 
# ripsTDA_diag <- ripsTDA$diagram[-1,]
# colnames(ripsTDA_diag) <- c('dimension', 'birth', 'death')
# plot_persist(ripsTDA_diag) + ggtitle('TDA')


#TDAstats
startTDAstats <- Sys.time()
ripsTDAstats <- calculate_homology(data, 
                                   dim = dimension, 
                                   threshold = xmax)
endTDAstats <-Sys.time()
print(c('TDAstats: ', endTDAstats-startTDAstats))
plot_persist(ripsTDAstats) + ggtitle('TDAstats')

