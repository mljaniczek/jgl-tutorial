library(igraph)
library(RColorBrewer)
library(scales)

# this is a funciton I made using code received from Kate Shutta and Yukun Li.
# Made this wrapper function to make the code within the tutorial document shorter.

plot_jgl <- function(
  thisGraph,
  multiplier = 15,
  vertex.size = 4,  
  vertex.label.dist = 2 , 
  vertex.label.cex = 0.5,
  rescale = F,
  vertex.color = NULL, 
  vertex.label.color= NULL, 
  main = NULL) {
  myLayout =  layout_in_circle(makePrettyGraphFromGraph(thisGraph))*0.8
  lab.locs <- radian.rescale(x=1:length(V(thisGraph)), direction=-1, start=0)
  plot(makePrettyGraphFromGraph(thisGraph, multiplier = multiplier), 
       vertex.size = vertex.size,  
       vertex.label.dist = vertex.label.dist , 
       vertex.label.cex = vertex.label.cex, 
       layout = myLayout,
       lab.locs = lab.locs,
       vertex.label.degree=lab.locs,
       vertex.color = vertex.color,
       vertex.label.color= vertex.label.color,
       rescale=F, 
       main = main)
}

# this is a function for making red-blue edges to match the edge weights(partial correlation)
# you may need to adjust the multiplier depending on how strong your edges are

makePrettyGraphFromGraph = function(thisGraph, multiplier = 15,redblue=T)
{
  if(redblue == T) my_palette <- brewer.pal(n = 8, name = "RdBu")
  if(redblue == F) my_palette <- brewer.pal(n = 8, name = "PRGn")
  my_palette <- my_palette[8:1]
  E(thisGraph)$width = abs(E(thisGraph)$weight)*multiplier
  E(thisGraph)$color = ifelse(sign(E(thisGraph)$weight)>0,my_palette[1],my_palette[8])
  return(thisGraph)
}


# this is code from this stack exchange: 
# https://stackoverflow.com/questions/23209802/placing-vertex-label-outside-a-circular-layout-in-igraph

radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

vBIC <- function(X, est_graph, thr=0.001){
  num <- length(X)
  BIC_acc <- 0.
  for(i in 1:num){
    
    data_num <- dim(X[[i]])[1]
    sample_cov <- cov(X[[i]], X[[i]])
    tr_sum <- sum(diag(sample_cov %*% est_graph$theta[[i]]))
    
    log_det <- determinant(est_graph$theta[[i]], logarithm = TRUE)$modulus[1][1]
    
    E <- sum(sum(abs(est_graph$theta[[i]]) >= thr))
    BIC_acc <- BIC_acc + (tr_sum - log_det) + (log(data_num)*E/data_num)
  }
  return(BIC_acc)
}
