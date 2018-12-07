library(TDA)
library(TDAstats)
library(tidyverse)

# example simulated 2D data
data <- cbind(runif(100), runif(100))
source('./Downloads/bettivector.R')

# step 0: preprocessing
# normalizing (like in cluster analysis, same problems when not)
data <- scale(data)

# step 1: convert data to diagram with TDA
# it results in a P by 3 matrix, 
# P is the number of points in the resulting persitence diagram
# the columns represents dimension, Birth, and Death of the features
# note: alpha_diag[1,3] == Inf

# with alpha complex
alpha_complex <- alphaComplexDiag(data)
alpha_diag <- alpha_complex$diagram

# with rips complex
rips_complex <- ripsDiag(data, maxdimension = 3, maxscale = 1)
rips_diag <- rips_complex$diagram

# step 2: convert diagram to betti curve
betticurve_data <- bettimatrix(alpha_complex)

# curvas de betti de grados 0 y 1
ggplot(betticurve_data, aes(x = xseq, y = value, color = order)) +
  geom_line() +
  facet_grid(rows = vars(order), scales = "free") + 
  ggtitle('titulo')

