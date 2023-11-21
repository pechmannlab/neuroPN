# Figure 1 neuroPN
# SP | 2022 | sebastian@pechmannlab.net

library(ggplot2)
library(reshape2)
library(ggdendro)
library(cowplot)

setwd("M2/neuro/")


# Figure 1A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


umap_data <- as.data.frame(read.table("data/processed/umap.txt", header=T))
colnames(umap_data) = c('umap1', 'umap2')

meta <- as.data.frame(read.table("data/processed/meta_types_cell.txt", header=T, sep='\t'))

umap_data$subtype <- meta$subtype
umap_data$subtype <- factor(umap_data$subtype, levels=sort(unique(umap_data$subtype)) )
list_clusters <- unique(umap_data$subtype)

cols_inh <- colorRampPalette(c("lightblue", "darkblue"))(length(grep("Inh", list_clusters)))
cols_exc <- colorRampPalette(c("yellow", "red"))(length(grep("Exc", list_clusters)))
cols_non <- colorRampPalette(c("grey80", "grey30"))(length(grep("non", list_clusters)))
colos <- c(cols_exc, cols_inh, cols_non)



svg(file = "figures/Figure1/A_umap.svg", height = 4, width = 4)

ggplot(umap_data) + 
  geom_point(aes(x=umap1, y=umap2, color=subtype), show.legend = F, size=0.01) +
  labs(x="UMAP1", y="UMAP2") + 
  scale_color_manual(values=colos) +
  theme_classic()+ 
  theme(
    text = element_text(size=18), 
    legend.title = element_blank()
  )

dev.off()






# Figure 1B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

exp_channels <- read.delim("data/processed/expression_channels.txt",  header=T, sep='\t')
exp_synapse <- read.delim("data/processed/expression_synapse.txt",  header=T, sep='\t')
exp_chap <- read.delim("data/processed/expression_chap.txt",  header=T, sep='\t')
exp_ub <- read.delim("data/processed/expression_ub.txt",  header=T, sep='\t')

meta.clusters <- read.delim("data/processed/meta_types_clusters.txt")




color_dendrogram <- function(exp_data, meta_data){
  
  # subroutine to compute cosine similarity
  cosinesimilarity <- function(a,b){
    d = sum(a*b)/sqrt(sum(a^2)*sum(b^2))
    return(d)
  }
  
  # data housekeeping
  row.names(exp_data) <- exp_data[,1]
  exp_data <- exp_data[,-1]
  exp_data <- as.matrix(exp_data)
  
  # cosine distance matrix
  cos.simil <- matrix(NA, ncol(exp_data), ncol(exp_data))
  labels.type = rep(NA, ncol(exp_data))
  for (i in 1:ncol(exp_data)){
    current_cluster <- colnames(exp_data)[i]
    current_idx <- which(meta_data$cluster==current_cluster)
    current_type = as.character(meta_data$type[current_idx])
    if (current_type != "Inh" & current_type != "Exc"){current_type <- 'non'}
    labels.type[i] <- current_type
    for (j in 1:ncol(exp_data)){
      s1 <- as.numeric(exp_data[,i])
      s2 <- as.numeric(exp_data[,j])
      sel <- is.na(s1) == F & is.na(s2) == F
      cos.simil[i,j] <- cosinesimilarity(s1[sel],s2[sel])
    }
  }

  # hierarchical clustering, adding labels, ...
  dist.cos <- as.dist(1-cos.simil)  # distance = 1 - similarity
  hc.cos <- hclust(dist.cos)
  hc.cos$labels <- labels.type
  cos.data <- dendro_data(hc.cos, type="rectangle")
  
  # DF for polygon coloring
  L <- length(cos.data$labels$label)
  id_counter <- 1
  ymax <- max(cos.data$segments$yend)
  y_vals <- c(0, ymax, ymax, 0)
  positions <- data.frame()
  datavals <- c()
  x1 <- 0.5
  for (i in 1:(L-1)){
    current <- cos.data$labels$label[i]
    prochain <- cos.data$labels$label[i+1]
    
    if (i < L-1){
      if (current == prochain){
        next
      }  
      if (current != prochain){
        x2 <- i+0.5
        x_vals <- c(x1, x1, x2, x2)
        poly <- data.frame(x=x_vals, y=y_vals, id=rep(factor(id_counter), 4))
        positions <- rbind(positions, poly)
        datavals <- c(datavals,as.character(current))
        id_counter <- id_counter + 1
        x1 <- x2
      }
    }
    if ( i == L-1){
      if (current == prochain){
        x2 <- L + 0.5
        x_vals <- c(x1, x1, x2, x2)
        poly <- data.frame(x=x_vals, y=y_vals, id=rep(factor(id_counter), 4))
        positions <- rbind(positions, poly)
        datavals <- c(datavals,as.character(current))
      }  
      if (current != prochain){
        x2 <- i+0.5
        x_vals <- c(x1, x1, x2, x2)
        poly <- data.frame(x=x_vals, y=y_vals, id=rep(factor(id_counter), 4))
        positions <- rbind(positions, poly)
        datavals <- c(datavals,as.character(current))
        id_counter <- id_counter + 1
        x1 <- x2
        
        x2 <- L + 0.5
        x_vals <- c(x1, x1, x2, x2)
        poly <- data.frame(x=x_vals, y=y_vals, id=rep(factor(id_counter), 4))
        positions <- rbind(positions, poly)
        datavals <- c(datavals,as.character(prochain))
      }
    }
  }
  
  values <- data.frame(id=unique(positions$id), value=datavals)
  datapoly <- merge(values, positions, by=c("id"))
  
  output <- list(cos.data, datapoly)
  return(output)
  
}


result_channels <- color_dendrogram(exp_channels, meta.clusters)
result_synapse <- color_dendrogram(exp_synapse, meta.clusters)
result_chap <- color_dendrogram(exp_chap, meta.clusters)
result_ub <- color_dendrogram(exp_ub, meta.clusters)

channels_data <- result_channels[[1]]
channels_poly <- result_channels[[2]]
synapse_data <- result_synapse[[1]]
synapse_poly <- result_synapse[[2]]
chap_data <- result_chap[[1]]
chap_poly <- result_chap[[2]]
ub_data <- result_ub[[1]]
ub_poly <- result_ub[[2]]



neurodendro <- function(exp_data, poly_data, title){

  y_max = max(exp_data$segments$yend)
  
  breaks <- c(0, 0.3, 0.6)
  if (y_max > 0.8){ breaks <- c(0, 0.5, 1)}
    
  p <- ggplot(segment(exp_data)) + 
    geom_polygon(data=poly_data, inherit.aes=F, aes(x=x, y=y, fill=value), alpha=0.5, show.legend=F) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend) ) + 
    scale_fill_manual(values=c("#2988E2", "yellow", "grey30")) + 
    coord_flip() + 
    scale_y_continuous(breaks=breaks, labels=breaks) + 
    theme_classic() + 
    labs(x="", y=bquote("d"["cos"]), title=title ) + 
    theme(
      axis.text.x = element_text(size=14),
      axis.title.x.bottom = element_text(size=18),
      axis.title.x.top = element_text(size=21),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
  return(p)
}


plot.channels <- neurodendro(channels_data, channels_poly, "Channels")
plot.synapse  <- neurodendro(synapse_data, synapse_poly, "Synapse")
plot.chap  <- neurodendro(chap_data, chap_poly, "Chaperones")
plot.ub  <- neurodendro(ub_data, ub_poly, "Ub ligases")


svg(file = "figures/Figure1/B_dendro.svg", height = 4, width = 6)

plot_grid(plot.channels, plot.synapse, plot.chap, plot.ub, labels ="", ncol = 4, align = 'h')

dev.off()

#ggdendrogram(cos.data, rotate=T, size=2) 


# Figure 1C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





diffexp <- as.data.frame(read.table("data/processed/classes_diffexp.txt", header=T, sep='\t'))
diffexp$class <- factor(diffexp$class, levels=c("Ch", "Syn", "Chap", "Ub"))         # keep order for plotting
diffexp$dir <- factor(diffexp$dir, levels=c("up", "down"))                          # keep order for plotting

svg(file = "figures/Figure1/C_diffexp.svg", height = 4, width = 2.2)

ggplot(diffexp) + 
  geom_col(aes(x=label, y=count, fill=dir), position = position_dodge2(), show.legend=F ) + 
  facet_grid(rows = vars(class)) + 
  scale_fill_manual(values=c("grey30", "grey80")) + 
  scale_y_continuous(breaks=c(0, 25, 50, 75), labels=c(0, 25, 50, 75)) + 
  labs(x="", y="% diff. expressed") + 
  theme_classic() + 
  theme(
    text = element_text(size=16),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=50, vjust = 0)
  )

dev.off()







# Figure 1D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# borrows variables from Fig 2!

var.chap <- as.matrix(read.table("data/processed/variability_chap.txt"))
exp.chap <- as.matrix(read.table("data/processed/expression_chap.txt")) 
all( colnames(var.chap) == colnames(exp.chap) )

meta.clusters <- read.delim("data/processed/meta_types_clusters.txt")



type <- c()
for (i in 1:ncol(var.chap)){
  current_cluster <- colnames(var.chap)[i]
  idx <-which(meta.clusters$cluster == current_cluster )
  current_type <- meta.clusters$type[idx] 
  current_type <- as.character(current_type)
  if (current_type != "Exc" & current_type != "Inh"){
    current_type <- "non"
  }
  type <- c(type, current_type )
}


cat.chap <- as.matrix(read.table("data/processed/abundance_bycategories_chap.txt"))
cat.ub <- as.matrix(read.table("data/processed/abundance_bycategories_ub.txt"))
cat.ch <- as.matrix(read.table("data/processed/abundance_bycategories_ch.txt"))
cat.syn <- as.matrix(read.table("data/processed/abundance_bycategories_syn.txt"))


df.all <- rbind(data.frame(total=colSums(cat.chap), type=type, cat=rep("Chap", length(type))), 
                data.frame(total=colSums(cat.ub), type=type, cat=rep("Ub", length(type))),  
                data.frame(total=as.numeric(cat.ch), type=type, cat=rep("Ch", length(type))), 
                data.frame(total=as.numeric(cat.syn), type=type, cat=rep("Syn", length(type))) )
df.all$cat <- factor(df.all$cat, levels=c("Ch", "Syn", "Chap", "Ub")) 


svg(file = "figures/Figure1/D_totalexp.svg", height = 4, width = 2.2)

ggplot(df.all) + 
  geom_boxplot(aes(x=type, y=total, fill=type), show.legend=F) + 
  facet_grid(rows = vars(cat), scales="free_y") +
  scale_fill_manual(values=c("yellow", "#2988E2", "grey50")) +
  labs(x="", y="Total expression") + 
  theme_classic() + 
  theme(
    text = element_text(size=14),
    axis.title.y = element_text(size=16), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size=16, angle=50, vjust = 0)
  )

dev.off()





