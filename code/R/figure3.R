# Figure 3 neuroPN
# SP | 2022 | sebastian@pechmannlab.net

library(ggplot2)
library(reshape2)
library(ggdendro)
library(cowplot)
library(igraph)
library(ggraph)


setwd("M2/neuro/")


# Figure 3A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


degree5 <- data.frame(read.table("data/processed/coexp_5_degree.txt", header=T))
degree5$class <- as.factor(degree5$class)


svg(file = "figures/Figure3/A_degree.svg", height = 3, width = 3)

ggplot(degree5, aes(x=degree, colour=class)) + 
  stat_ecdf(size=1.5) + 
  labs(x="Degree", y="Fraction") + 
  scale_x_continuous(limits=c(0, 2200)) + 
  scale_color_manual(values=c("grey50", "orange", "purple")) +
  coord_flip() + 
  theme_classic() + 
  theme(
    text = element_text(size=18),
    legend.position = c(0.9, 0.25),
    legend.title = element_blank()
  )

#ggplot(degree5) + 
#  geom_boxplot(aes(x=class, y=degree, fill=class), show.legend=F) +
#  labs(x="", y="Degree") +
#  scale_fill_manual(values=c("grey50", "orange", "purple")) +
#  theme_classic() + 
#  theme(
#    text = element_text(size=18), 
#    #axis.text.x = element_text(angle=60, vjust=1),
#    axis.line.x = element_blank(),
#    axis.ticks.x = element_blank()
#  )

dev.off()




# Figure 3B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

degree5 <- data.frame(read.table("data/processed/coexp_5_degree.txt", header=T))
deg_exc_inh <- data.frame(read.table("data/processed/coexp_exc_inh_degree.txt", header=T))
deg_inh_exc <- data.frame(read.table("data/processed/coexp_inh_exc_degree.txt", header=T))
deg_neu_non <- data.frame(read.table("data/processed/coexp_neu_non_degree.txt", header=T))



chap.counts <- data.frame(cat=c("Total", "Neuron vs non", "Exc vs Inh", "Inh vs Exc"),
                        counts=c(sum(degree5$class=="Chap"), 
                                 sum(deg_neu_non[deg_neu_non$class=="Chap",2] > 0),
                                 sum(deg_exc_inh[deg_exc_inh$class=="Chap",2] > 0),
                                 sum(deg_inh_exc[deg_inh_exc$class=="Chap",2] > 0))
)


ub.counts <- data.frame(cat=c("Total", "Neuron vs non", "Exc vs Inh", "Inh vs Exc"),
                        counts=c(sum(degree5$class=="Ub"), 
                                 sum(deg_neu_non[deg_neu_non$class=="Ub",2] > 0),
                                 sum(deg_exc_inh[deg_exc_inh$class=="Ub",2] > 0),
                                 sum(deg_inh_exc[deg_inh_exc$class=="Ub",2] > 0))
                        )
                             

plot.chap <- ggplot(chap.counts) + 
  geom_col(aes(x=cat, y=counts), fill="orange") +
  labs(x="", y="Count") +
  coord_flip() +
  scale_y_continuous( breaks=c(0, 100, 200), labels=c(0, 100, 200) ) +
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

plot.ub <- ggplot(ub.counts) + 
  geom_col(aes(x=cat, y=counts), fill="purple") +
  labs(x="", y="Count") +
  coord_flip() +
  scale_y_continuous( breaks=c(0, 150, 300), labels=c(0, 150, 300) ) +
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )



svg(file = "figures/Figure3/B_counts.svg", height = 3, width = 3)

plot_grid(plot.chap, plot.ub, labels ="", ncol = 1, align = 'v')

dev.off()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 3C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vertices.neu_non <- as.data.frame(read.table("data/processed/clusters/neu_non_annotation.txt", header=T))

ind <- 0
df <- data.frame(Gene=NULL, cat=NULL, cluster=NULL, degree=NULL, PNsyst=NULL, pos=NULL)
for (i in 1:length(levels(vertices.neu_non$cluster))){
  current_vert <- vertices.neu_non[vertices.neu_non$cluster==levels(vertices.neu_non$cluster)[i], ]
  current_vert <- current_vert[current_vert$cat!="syn", ]
  current_vert <- current_vert[sort(current_vert$degree, index.return=T, decreasing=T)$ix,]
  current_vert <- current_vert[current_vert$degree > 20,]
  current_vert <- rbind(current_vert, data.frame(Gene="", cat="X", cluster=levels(vertices.neu_non$cluster)[i], degree=-45, PNsyst="None"))

    if (nrow(current_vert) > 0) {
    current_vert$pos <- c(1:nrow(current_vert)) + ind
    #print(current_vert)
    df <- rbind(df, current_vert)
    ind <- max(current_vert$pos) + 1
    }
  df
}

tmp <- df$Gene == ""
i <- 0
keep <- rep(T, length(tmp))
while (i >= 0){
  n = length(tmp) - i
  if (tmp[n] == T){
    print(i)
    keep[n] = F
    i <- i + 1
  }
  else { i <- - 100}
}
df <- df[keep,]

svg("figures/Figure3/C_degree.svg", height=6, width=2)

ggplot(df) + 
  geom_col(aes(x=pos, y=degree, fill=cat), show.legend=F) + 
  geom_text(aes(x=pos, y=-0.5, label=Gene), hjust=1, size=3) + 
  scale_fill_manual(values=c("orange", "purple", "white")) +
  labs(x="", y="Degree") + 
  coord_flip() + 
  scale_x_reverse() + 
  theme_classic() +
  theme(
    text = element_text(size=16),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 3D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


network <- data.frame(read.table("data/processed/clusters/neu_non_network.txt", header=T))

network$alpha <- ifelse(network$cluster == "noclust", 0.1, 1)
network$width <- ifelse(network$cluster == "noclust", 0.01, 0.06)

vertices <- as.data.frame(read.table("data/processed/clusters/neu_non_annotation.txt", header=T))

graph <- graph_from_data_frame(network, vertices=as.vector(vertices$Gene))
net2 <- network[network$cluster!="nocluster",]
graph2 <- graph_from_data_frame(net2, vertices=as.vector(vertices$Gene))

vertices$degree <- degree(graph)
vertices$degree2 <- degree(graph2)
vertices$size <- ifelse(vertices$cat != "syn", vertices$degree, 1)

c1 <- colorRampPalette(c("#29ABE2", "darkblue"))(5)
lay <- create_layout(graph, layout='igraph', algorithm='kk')  # create random layout and then change coords manually
lay <- as.matrix(lay[,c(1,2)])


svg("figures/Figure3/D_network.svg", height=5, width=5)

ggraph(graph2, layout=lay) +
  geom_edge_link(aes(colour=as.factor(net2$cluster), edge_alpha=net2$alpha), width=0.03 ) + 
  scale_edge_color_manual(values=c(c1, 'grey90')) + 
  geom_node_point(aes(color=factor(vertices$cat), size=vertices$size) ) +
  scale_color_manual(values=c( "orange", "grey50", "purple")) +  
  #geom_node_text(aes(label=vertices$Gene)) + 
  theme_void() + 
  theme(
    legend.position='none'
  )
dev.off()
  




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 3E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

deg <- data.frame(read.table("data/processed/clusters/degree.txt", header=T))
deg$index <- factor(c("Exc vs Inh", "Inh vs Exc", "Neuron vs non"), levels=c("Neuron vs non", "Exc vs Inh", "Inh vs Exc"))


svg("figures/Figure3/E_reldegree.svg", height=4, width=2.5)

ggplot(deg) + 
  geom_col(aes(x=index, y=log2(OR), fill=index), show.legend=F ) +  
  labs(x="", y="Rel. average degree") + 
  scale_fill_manual(values=c(  )) + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(angle=40, hjust=1),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )

dev.off()
 

syst <- data.frame(read.table("data/processed/clusters/degree_system.txt", header=T))
syst <- syst[which(rowSums(syst[,c(3,5,7)] < 0.05) > 0),]       # not show those without any significance

syst_sig <- syst[,c(1,3,5,7)]
syst_dat <- syst[,c(1,2,4,6)]
syst_m <- melt(syst_dat, id="System")
syst_p <- melt(syst_sig, id="System")
syst_m$pval <- syst_p$value
syst_m$alpha <- ifelse(syst_m$pval < 0.05, 1, 0.5)
syst_m$variable <- factor(syst_m$variable, levels=c("enrich_Nn", "enrich_EI", "enrich_IE"))
syst_m$System <- factor(syst_m$System, levels=c("CCT/TRiC_system", "HSP70_system", "HSP90_system", "RING_Single", "RING_Complex"))



svg("figures/Figure3/E_catdegree.svg", height=4, width=3)

ggplot(syst_m) + 
  geom_col(aes(x=System, y=log2(value), fill=variable, alpha=alpha) , position = position_dodge2(), show.legend=F) + 
  labs(x="", y="Rel. average degree") + 
  scale_fill_manual(values=c( "#222222", "#777777", "#BBBBBB")) + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(angle=40, hjust=1),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 3F ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


network_exc <- data.frame(read.table("data/processed/clusters/exc_inh_network.txt", header=T))
vertices_exc <- as.data.frame(read.table("data/processed/clusters/exc_inh_annotation.txt", header=T))


network_exc$alpha <- ifelse(network_exc$cluster == "cluster0", 1, 0.2)
network_exc$width <- ifelse(network_exc$cluster == "cluster0", 0.1, 0.02)

graph_exc <- graph_from_data_frame(network_exc, vertices=as.vector(vertices_exc$Gene))
net2_exc <- network_exc[network_exc$cluster!="nocluster",]
graph2_exc <- graph_from_data_frame(net2_exc, vertices=as.vector(vertices_exc$Gene))

vertices_exc$degree <- degree(graph_exc)
vertices_exc$degree2 <- degree(graph2_exc)
vertices_exc$size <- ifelse(vertices_exc$cat != "syn", vertices_exc$degree, 1)

c1 <- colorRampPalette(c("#29ABE2", "darkblue"))(5)
lay <- create_layout(graph_exc, layout='igraph', algorithm='kk')  # create random layout and then change coords manually
lay <- as.matrix(lay[,c(1,2)])


svg("figures/Figure3/F_network_all.svg", height=4, width=4)

ggraph(graph2_exc, layout=lay) +
  geom_edge_link(aes(colour=as.factor(net2_exc$cluster), edge_alpha=net2_exc$alpha), width=0.08 ) + 
  scale_edge_color_manual(values=c(c1, 'grey90')) + 
  geom_node_point( aes(color=factor(vertices_exc$cat)), size=0.7 ) +
  scale_color_manual(values=c( "orange", "grey50", "purple")) +  
  #geom_node_text(aes(label=vertices$Gene)) + 
  theme_void() + 
  theme(
    legend.position='none'
  )
  
dev.off()









network_exc <- data.frame(read.table("data/processed/clusters/exc_inh_network.txt", header=T))
vertices_exc <- as.data.frame(read.table("data/processed/clusters/exc_inh_annotation.txt", header=T))


n0 <- network_exc[network_exc$cluster == "cluster0",]
v0 <- vertices_exc[vertices_exc$cluster == "cluster0",]
graph <- graph_from_data_frame(n0, vertices=as.vector(v0$Gene))

lay <- create_layout(graph, layout='igraph', algorithm='kk')  # create random layout and then change coords manually
lay <- as.matrix(lay[,c(1,2)])


svg("figures/Figure3/F_network_cluster.svg", height=4, width=4)

ggraph(graph, layout=lay) +
  geom_edge_link(color="#29ABE2", width=0.1 ) + 
  #scale_edge_color_manual(values=c(c1, 'grey90')) + 
  geom_node_point(aes(color=factor(v0$cat), size=v0$degree ) ) +
  scale_color_manual(values=c( "orange", "grey50", "purple")) +  
  #geom_node_text(aes(label=v0$Gene)) + 
  theme_void() + 
  theme(
    legend.position='none'
  )

dev.off()









v_hsp90 <- vertices_exc[vertices_exc$PNsyst=="HSP90_system",]

list_nodes <- as.character(v_hsp90$Gene)
for (i in 1:nrow(v_hsp90)){
  current_hsp <- as.character(v_hsp90$Gene[i])
  current_neighbors <- neighbors(graph, current_hsp)
  list_nodes <- c(list_nodes, names(current_neighbors))
}

graph_hsp90 <- induced_subgraph(graph, list_nodes, delete.vertices=F)
list_v <- names(V(graph_hsp90))

list_idx <- c()
for (i in 1:length(list_v)){
  current_node <- list_v[i]
  current_v <- which(vertices_exc$Gene==current_node)
  list_idx <- c(list_idx, current_v)
}
vertices_hsp90 <- vertices_exc[list_idx,]
vertices_hsp90$size <- ifelse(vertices_hsp90$PNsyst == "HSP90_system", 2, 2)
vertices_hsp90$label <- ifelse(vertices_hsp90$PNsyst == "HSP90_system", as.character(vertices_hsp90$Gene), "")

lay <- create_layout(graph_hsp90, layout='igraph', algorithm='fr')  # create random layout and then change coords manually
lay <- as.matrix(lay[,c(1,2)])


svg("figures/Figure3/F_network_hsp90.svg", height=4, width=4)

ggraph(graph_hsp90, layout=lay) +
  geom_edge_link(width=0.1 ) + 
  #scale_edge_color_manual(values=c(c1, 'grey90')) + 
  geom_node_point( aes(color=vertices_hsp90$cat, size=vertices_hsp90$degree) ) +
  scale_color_manual(values=c( "orange", "grey50", "purple")) +  
  geom_node_text(aes(label=vertices_hsp90$label,  x=x*1.15, y=y*1.15), size=5, alpha=1, repel=T) + 
  theme_void() + 
  theme(
    legend.position='none'
  )

dev.off()



