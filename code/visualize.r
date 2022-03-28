library(igraph)
library(RColorBrewer)

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

library(scales)

radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }

lab.locs <- radian.rescale(x=1:length(V(thisGraph)), direction=-1, start=0)

# multiplying the layout * 0.8 shrinks everything slightly
# layout is just coordinates in x,y plane so you can 
# multiply it by any scalar you want

myLayout = layout_in_circle(makePrettyGraphFromGraph(thisGraph))*0.8

# here you can assign the color of the vertex however you like
V(thisGraph)$color = ifelse(V(thisGraph)$name %in% discOnly, "blue", ifelse(V(thisGraph)$name %in% val, "red","darkgreen"))

plot(makePrettyGraphFromGraph(thisGraph), 
     vertex.size = 4,  
     vertex.label.dist = 2 , 
     vertex.label.cex = 0.5, 
     vertex.label.color=V(thisGraph)$color,  
     layout = myLayout, vertex.label.degree=lab.locs,rescale=F)

#######
# myGraph = graph_from_adjacency_matrix(-cov2cor(myMatrix), 
#                                       weighted=T,mode="undirected",diag=F)
# plot(myGraph)

thisGraph = graph_from_adjacency_matrix(-cov2cor(cov.mat), 
                                      weighted=T,mode="undirected",diag=F)

thisGraph = graph_from_adjacency_matrix(-cov2cor(inv.mat), 
                                      weighted=T,mode="undirected",diag=F)

thisGraph = graph_from_adjacency_matrix(-cov2cor(cov.mat.X), 
                                      weighted=T,mode="undirected",diag=F)

thisGraph = graph_from_adjacency_matrix(-cov2cor(inv.mat.X), 
                                      weighted=T,mode="undirected",diag=F)

lab.locs <- radian.rescale(x=1:length(V(thisGraph)), direction=-1, start=0)

myLayout = layout_in_circle(makePrettyGraphFromGraph(thisGraph))*0.8

plot(makePrettyGraphFromGraph(thisGraph, multiplier = 2), 
     vertex.size = 4,  
     vertex.label.dist = 2 , 
     vertex.label.cex = 0.5, 
     vertex.label.color=V(thisGraph)$color,  
     layout = myLayout, vertex.label.degree=lab.locs,rescale=F)
