######################################## SETTINGS ######################################## 

## Package load
library(igraph)
library(intergraph)
library(statnet)
library(sna)
library(rms)
library(kableExtra)
library(knitr)
library(moments)

## Directory set
setwd('./Validation/data/old')

######################################## DATA CODING #########################################

Sampson_paj <- read.paj("Sampson.paj")
Sampson_T4 <- Sampson_paj$networks$Sampson_T4
Sampson_net <- Sampson_paj$networks$Sampson

## Graph and Network creation ###
Sampson <- asIgraph(Sampson_net)
Sampson_graph <- asIgraph(Sampson_T4)

## Main Nodelist creation ###
Sampson_nodes <- as.data.frame(1:25)
Sampson_nodes[,2] <- V(Sampson)$vertex.names
Sampson_nodes[,3] <- Sampson_paj$partitions$Sampson_cloisterville
colnames(Sampson_nodes) <- c("Id", "Label", "Cloisterville")

## Nodelist creation from T4 Network ###
Sampson_nodes_T4 <- as.data.frame(1:18)
Sampson_nodes_T4[,1] <- V(Sampson_graph)$vertex.names
Sampson_nodes_T4[,2] <- Sampson_paj$partitions$Sampson_cloisterville_T4
Sampson_nodes_T4[,3] <- Sampson_paj$partitions$Sampson_factions_T4
colnames(Sampson_nodes_T4) <- c("Label", "Cloisterville", "Factions")

# Adding the T4 data to the main network
Sampson_nodes <- merge(Sampson_nodes, 
                       Sampson_nodes_T4, 
                       by = c("Label", "Cloisterville"), 
                       all = TRUE)
Sampson_nodes <- Sampson_nodes[c(3,1,2,4)]

## Edgelist creation ###
# With the elimination of repeated links
Sampson_edges <- as.data.frame(as.edgelist(Sampson_net, attrname = "Sampson"))
colnames(Sampson_edges) <- c("Source", "Target", "Relation")
Sampson_edges_pos <- Sampson_edges[which(Sampson_edges$Relation > 0), 1:2]
Sampson_edges_pos <- unique(Sampson_edges_pos, by = "Target")
Sampson_edges_neg <- Sampson_edges[which(Sampson_edges$Relation < 0), 1:2]
Sampson_edges_neg <- unique(Sampson_edges_neg, by = "Target")
Sampson_edges$Relation <- ifelse(Sampson_edges$Relation > 0, "Positive", "Negative")
Sampson_edges <- unique(Sampson_edges, by = "Target")

## Data export ###
write.csv2(Sampson_nodes, file = "Sampson_nodes.csv", row.names=FALSE)
write.csv2(Sampson_edges, file = "Sampson_edges.csv", row.names=FALSE)
write.csv2(Sampson_edges_neg, file = "Sampson_edges_neg.csv", row.names=FALSE)
write.csv2(Sampson_edges_pos, file = "Sampson_edges_pos.csv", row.names=FALSE)

######################################## EXCERCICE######################################## 
## Data Import ####
setwd('./Validation/data')

Snodes <- data.frame(read.csv2("Sampson_nodes.csv", header=T, as.is=T))
Sedges <- data.frame(read.csv2("Sampson_edges.csv", header=T, as.is=T))
Sedges_pos <- data.frame(read.csv2("Sampson_edges_pos.csv", header=T, as.is=T))
Sedges_neg <- data.frame(read.csv2("Sampson_edges_neg.csv", header=T, as.is=T))

## Graph creation ####
Sall <- graph.data.frame(Sedges, Snodes, directed=TRUE)
Sneg <- graph.data.frame(Sedges_neg, Snodes, directed=TRUE)
Spos <- graph.data.frame(Sedges_pos, Snodes, directed=TRUE)

## Graph plot ####
Sall_layout <- layout.fruchterman.reingold(Sall)
plot(Sall,
     edge.arrow.size = .1,
     rescale = TRUE,
     vertex.color = V(Sall)$Factions,
     vertex.shape = ifelse(V(Sall)$Cloisterville == "1", "square","circle"),  
     edge.color= ifelse(E(Sall)$Relation == "Positive","blue","red"))

plot(Spos,
     edge.arrow.size = .1,
     rescale = TRUE,
     vertex.color = V(Sall)$Factions,
     vertex.shape = ifelse(V(Sall)$Cloisterville == "1", "square","circle"),  
     edge.color = "blue")

plot(Sneg,
     edge.arrow.size = .1,
     rescale = TRUE,
     vertex.color = V(Sall)$Factions,
     vertex.shape = ifelse(V(Sall)$Cloisterville == "1", "square","circle"),  
     edge.color = "red")

## CENTRALITY CALCULATION ####
## In-Degree Centrality #####
# In-degre centrality (pos. influence)
Spos_deg_in <- centr_degree(Spos, mode = "in")
Snodes[,5:8] <- Spos_deg_in
# In-degre centrality (neg. influence)
Sneg_deg_in <- centr_degree(Sneg, mode = "in")

## Betweeness Centrality ####
Spos_bet <- igraph::betweenness(Spos, directed=TRUE)
Sneg_bet <- igraph::betweenness(Sneg, directed=TRUE)

## Closeness Centrality ####
Spos_close_in <- igraph::closeness(Spos, mode='in')
Sneg_close_in <- igraph::closeness(Sneg, mode='in')

## Eigenvector Centrality ####
Sneg_un <- graph.data.frame(Sedges_neg, Snodes, directed=FALSE)
Spos_un <- graph.data.frame(Sedges_pos, Snodes, directed=FALSE)

Spos_eigen <- eigen_centrality(Spos_un, directed = FALSE, scale = TRUE, weights = NULL)
Sneg_eigen <- eigen_centrality(Sneg_un, directed = FALSE, scale = TRUE, weights = NULL)

# Centrality comparison w/ other centrality measures
Sall_cen <- data.frame(Spos_deg_in$res, Spos_close_in, Spos_bet, Spos_eigen$vector,
                       Sneg_deg_in$res, Sneg_close_in, Sneg_bet, Sneg_eigen$vector)
Sall_cen_pos <- data.frame(Spos_deg_in$res, Spos_close_in, Spos_bet, Spos_eigen$vector)
Sall_cen_neg <- data.frame(Sneg_deg_in$res, Sneg_close_in, Sneg_bet, Sneg_eigen$vector)
car::scatterplotMatrix(~Spos_deg_in$res+Spos_close_in+Spos_bet+Spos_eigen$vector,
      data = Sall_cen_pos)

# Centrality comparison w/ Factions
Sall_fac <- data.frame(V(Sall)$Factions,
                       Spos_deg_in$res, Spos_close_in, Spos_bet,
                       Sneg_deg_in$res, Sneg_close_in, Sneg_bet)
Sall_fac[is.na(Sall_fac)] <- 0

summary(Sall_fac) # I transcribed the summary data into a new csv table:
Spos_cen <- read.csv2("degree_pos.csv")
Sneg_cen <- read.csv2("degree_neg.csv")

## NETWORK METRICS ####
## Density #####
graph.density(Spos)
graph.density(Sneg)

## Av. Path Length (Connectivity) ####
# Distance moyenne entre les noeuds
average.path.length(Spos)

# Diameter
diameter(Spos, weights = NA)

# Eccentricity
Spos_ecce <- eccentricity(Spos)
Sneg_ecce <- eccentricity(Sneg)

# Radius
radius(Spos)

## Transitivity #####
# Global
transitivity(Spos)
transitivity(Sneg)
# Local
Spos_tran <- transitivity(Spos, type = "local")
Sneg_tran <- transitivity(Sneg, type = "local")

transitivity(Spos, type = "local")
transitivity(Spos, type = "localaverage")
transitivity(Sneg)

# Transitivity Comparison
Sall_tran <- data.frame(V(Sall)$Factions, Spos_tran, Sneg_tran)
Sall_tran[is.na(Sall_tran)] <- 0

summary(Sall_tran) # I transcribed the summary data into a new csv table:
Sall_tran <- read.csv2("transitivity.csv")

# Measures of skewness for each faction
skewness(Sall_fac[Sall_fac$V.Sall..Factions %in% 0,])
skewness(Sall_fac[Sall_fac$V.Sall..Factions %in% 1,])
skewness(Sall_fac[Sall_fac$V.Sall..Factions %in% 2,])
skewness(Sall_fac[Sall_fac$V.Sall..Factions %in% 3,])
skewness(Sall_fac[Sall_fac$V.Sall..Factions %in% 4,])
