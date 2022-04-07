library(ggraph)

test_graph <- graph_from_adjacency_matrix(
  -cov2cor(inv_covar_matrices[[1]]),
  weighted = T,
  mode = "undirected",
  diag = FALSE
)

ggraph(test_graph) +
  geom_edge_link() + 
  geom_node_point()

library(GGally)
library(network)
library(sna)
library(intergraph)

color_opts <- RColorBrewer::brewer.pal(3, "Dark2")
color_list <- data.frame(class = met_group) %>%
  mutate(color = case_when(
    class == "Amino Acids" ~ color_opts[1],
    class == "Acyl carnitines" ~ color_opts[2],
    class == "Other" ~ color_opts[3],
    TRUE ~ NA_character_
  )
  )

multiplier = 8
my_palette <- brewer.pal(n = 8, name = "RdBu")
my_palette <- my_palette[8:1]
E(test_graph)$width = abs(E(test_graph)$weight)*multiplier
E(test_graph)$color = ifelse(sign(E(test_graph)$weight)>0,my_palette[1],my_palette[8])


ggnet2(test_graph,
       mode = "circle",
       color = color_list$color,
       edge.color = E(test_graph)$color,
       edge.size = abs(E(test_graph)$weight)*2,
       #label = TRUE,
       label.size = 3,
       #vjust = -1,
       size = 2)
