# Figure 2 neuroPN
# SP | 2022 | sebastian@pechmannlab.net

library(ggplot2)
library(reshape2)
library(ggdendro)
library(cowplot)
library(viridis)





# Figure 2A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


c1 <- colorRampPalette(c("#29ABE2", "white"))(7)
c2 <- colorRampPalette(c("white", "gold"))( 7 )
col <- c(c1[c(1,3,6)], c2[c(1,2,5,7)])

#col <- c("#29ABE2", "#70C7EB", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFE455", "#FFD700")
#col <- c("#29ABE2", "#70C7EB", "#222222", "#222222", "#222222", "#FFE455", "#FFD700")
#col <- c("#29ABE2", "#70C7EB", "#FFFFFF", "#FFE455", "#FFD700")

col <- c("#29ABE2", "#333333", "#FFD700", "#FFFFFF")

exp.chap <- read.table("data/processed/heatmap_chap.txt", header=T, sep='\t')
colnames(exp.chap)[1] <- "Gene"
exp.chap$Gene <- factor(exp.chap$Gene, levels=exp.chap$Gene)
exp.chap <- as.data.frame(exp.chap)
data_chap <- melt(exp.chap, id="Gene")
data_chap$log <- log2(data_chap$value) #log2
data_chap$log[data_chap$log == -Inf] <- 0
data_chap$log[is.na(data_chap$log)] <- 5000
#data_chap$breaks <- cut(as.numeric(data_chap$log), c(-1000,-2, -1, -0.5, 0.5, 1, 2, 1000), right=F)
data_chap$breaks <- cut(as.numeric(data_chap$log), c(-1000, -1,  1, 1000, 10000), right=F)


exp.ub <- read.table("data/processed/heatmap_ub.txt", header=T, sep='\t')
colnames(exp.ub)[1] <- "Gene"
exp.ub$Gene <- factor(exp.ub$Gene, levels=exp.ub$Gene)
exp.ub <- as.data.frame(exp.ub)
data_ub <- melt(exp.ub, id="Gene")
data_ub$log <- log2(data_ub$value)
data_ub$log[data_ub$log == -Inf] <- 0
data_ub$log[is.na(data_ub$log)] <- 5000
#data_ub$breaks <- cut(as.numeric(data_ub$log), c(-1000,-2, -1, -0.5, 0.5, 1, 2, 1000), right=F)
data_ub$breaks <- cut(as.numeric(data_ub$log), c(-1000, -1,  1, 1000, 10000), right=F)


plot.chap <- ggplot(data_chap, aes(x = variable, y = Gene)) + 
  geom_tile(aes(fill=breaks),color="white", size=0.01) +
  scale_fill_manual(values=col) +
  labs(x=NULL, y="Chaperones") + 
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size=16),
    axis.ticks=element_blank() ) 


plot.ub <- ggplot(data_ub, aes(x = variable, y = Gene)) + 
  geom_tile(aes(fill=breaks),color="white", size=0.01) +
  scale_fill_manual(values=col) +
  labs(x=NULL, y="Ub ligases") + 
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size=16),
    axis.ticks=element_blank() ) 


#svg(file = "figures/Figure2/A_heat.svg", height = 5, width = 5)

png(filename="figures/Figure2/A_heat.png", width=5, height=5, unit="in", res=1200)

plot_grid(plot.chap, plot.ub, labels ="", ncol = 1, align = 'v')

dev.off()



# Figure 2B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


l2fc.chap <- as.data.frame(read.table("data/processed/log2fc_sig_chap.txt", header=T, sep='\t'))
top10.chap <- l2fc.chap[ sort(rowSums(l2fc.chap[,c(2:7)]), index.return=T)$ix[(nrow(l2fc.chap)-9):nrow(l2fc.chap)], ]
top10.chap <- melt(top10.chap, id="X")
top10.chap$cat <- c(rep("Exc", 20), rep("Inh", 20), rep("non", 20))   # need to check that still matchings labels after changes

top4.chap <- l2fc.chap[ sort(rowSums(l2fc.chap[,c(2:7)]), index.return=T)$ix[(nrow(l2fc.chap)-3):nrow(l2fc.chap)], ]
top4.chap <- melt(top4.chap, id="X")
top4.chap$cat <- c(rep("Exc", 8), rep("Inh", 8), rep("non", 8))   # need to check that still matchings labels after changes

examples.chap <- l2fc.chap[c(which(l2fc.chap$X=="DNAJB13"), which(l2fc.chap$X=="FKBP5"), which(l2fc.chap$X=="HSPB9"), which(l2fc.chap$X=="ZMYND10")),]
examples.chap <- melt(examples.chap, id="X")
examples.chap$cat <- c(rep("Exc", 8), rep("Inh", 8), rep("non", 8))   # need to check that still matchings labels after changes


* Hsp40 ER
DNAJC5G Hsp40 ER
DNAJC5B Hsp40 ER
*FKBP5 coHsp90
*HSPB9 sHsp
ODF1 sHsp
*ZMYND10 coHsp90

svg(file = "figures/Figure2/B_2fold_chap.svg", height = 3.5, width = 3)

ggplot(examples.chap) + 
  geom_col(aes(x=cat, y=value, fill=variable), position=position_dodge2(),  show.legend=F) + 
  facet_grid(rows = vars(X)) +  # , scales='free_y'
  scale_fill_manual(values=c("gold", "blue", "gold", "blue", "gold", "blue")) + 
  scale_y_continuous(labels=c(0, 50, 100), breaks=c(0, 50, 100)) + 
  labs(x="", y="% with > 2-fold change") + 
  theme_classic() + 
  theme(
    axis.title.y = element_text(size=20),
    axis.text.y = element_text(size=12), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.text.x = element_text(size=16)
    )

dev.off()



l2fc.ub <- as.data.frame(read.table("data/processed/log2fc_sig_ub.txt", header=T, sep='\t'))
top10.ub <- l2fc.ub[ sort(rowSums(l2fc.ub[,c(2:7)]), index.return=T)$ix[(nrow(l2fc.ub)-9):nrow(l2fc.ub)], ]
top10.ub <- melt(top10.ub, id="X")
top10.ub$cat <- c(rep("Exc", 20), rep("Inh", 20), rep("non", 20))   # need to check that still matchings labels after changes

top4.ub <- l2fc.ub[ sort(rowSums(l2fc.ub[,c(2:7)]), index.return=T)$ix[(nrow(l2fc.ub)-3):nrow(l2fc.ub)], ]
top4.ub <- melt(top4.ub, id="X")
top4.ub$cat <- c(rep("Exc", 8), rep("Inh", 8), rep("non", 8))   # need to check that still matchings labels after changes

examples.ub <- l2fc.ub[c(which(l2fc.ub$X=="TRIM22"), which(l2fc.ub$X=="RNF135"), which(l2fc.ub$X=="FBXO40"), which(l2fc.ub$X=="KCTD11")),]
examples.ub <- melt(examples.ub, id="X")
examples.ub$cat <- c(rep("Exc", 8), rep("Inh", 8), rep("non", 8))   # need to check that still matchings labels after changes


svg(file = "figures/Figure2/B_2fold_ub.svg", height = 3.5, width = 3)

ggplot(examples.ub) + 
  geom_col(aes(x=cat, y=value, fill=variable), position=position_dodge2(),  show.legend=F) + 
  facet_grid(rows = vars(X)) +  # , scales='free_y'
  scale_fill_manual(values=c("gold", "blue", "gold", "blue", "gold", "blue")) + 
  scale_y_continuous(labels=c(0, 40, 80), breaks=c(0, 40, 80)) + 
  labs(x="", y="% with > 2-fold change") + 
  theme_classic() + 
  theme(
    axis.title.y = element_text(size=20),
    axis.text.y = element_text(size=12), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.text.x = element_text(size=16)
  )

dev.off()




# Figure 2C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# variability

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

sel.chap <- rowSums( var.chap == 0 ) < 130 # exclude genes that are never variable in these clusters
#boxplot( rowMeans(var.chap[sel.chap, type=='Exc']), rowMeans(var.chap[sel.chap, type=='Inh']), rowMeans(var.chap[sel.chap, type=='non']), ylim=c(0, 0.4), col=2)

df.chap <- data.frame( Exc=rowMeans(var.chap[sel.chap, type=='Exc']), Inh=rowMeans(var.chap[sel.chap, type=='Inh']), non=rowMeans(var.chap[sel.chap, type=='non'] ) )
df.chap <- melt(df.chap)


svg(file = "figures/Figure2/C_var_chap.svg", height = 4, width = 3)

ggplot(df.chap) + 
  geom_boxplot(aes(x=variable, y=value, fill=variable), show.legend=F) + 
  labs(x="", y="rel% variable") + 
  scale_fill_manual(values=c("yellow", "#2988E2", "grey50")) + 
  scale_y_continuous(limits=c(0, 0.3), breaks=c(0, 0.1, 0.2, 0.3), labels=c(0, 10, 20, 30)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )

dev.off()


wilcox.test( df.chap[df.chap$variable=="Exc", 2], df.chap[df.chap$variable=="non", 2] )




var.ub <- as.matrix(read.table("data/processed/variability_ub.txt"))
sel.ub <- rowSums( var.ub == 0 ) < 130 # exclude genes that are never variable in these clusters
df.ub <- data.frame( Exc=rowMeans(var.ub[sel.ub, type=='Exc']), Inh=rowMeans(var.ub[sel.ub, type=='Inh']), non=rowMeans(var.ub[sel.ub, type=='non'] ) )
df.ub <- melt(df.ub)


svg(file = "figures/Figure2/C_var_ub.svg", height = 4, width = 3)

ggplot(df.ub) + 
  geom_boxplot(aes(x=variable, y=value, fill=variable), show.legend=F) + 
  labs(x="", y="rel% variable") + 
  scale_fill_manual(values=c("yellow", "#2988E2", "grey50")) + 
  scale_y_continuous(limits=c(0, 0.3), breaks=c(0, 0.1, 0.2, 0.3), labels=c(0, 10, 20, 30)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )

dev.off()





df2 <- data.frame( exc=rowMeans(exp.chap[sel, type=='Exc']), inh=rowMeans(exp.chap[sel, type=='Inh']), non=rowMeans(exp.chap[sel, type=='non'] ) )



# Figure 2D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# totals, multiple cats


cat.chap <- as.matrix(read.table("data/processed/abundance_bycategories_chap.txt"))
cat.chap <- cat.chap[c("HSP90_system", "HSP70_system", "CCT/TRiC_system", "sHSP_system"),]
row.names(cat.chap) <- c("HSP90", "HSP70", "TRiC", "sHSP")
all( colnames(var.chap) == colnames(cat.chap) )       # check that cols are corresponding
cat.chap <- as.data.frame(t(cat.chap))
cat.chap$type <- type
cat.chap.m <- melt(cat.chap, id="type")
           

svg(file = "figures/Figure2/D_total_chap.svg", height = 3.5, width = 4)

ggplot(cat.chap.m) + 
  geom_boxplot(aes(x=type, y=value, fill=type), show.legend = F ) + 
  labs(x="", y="Cumulative expression") + 
  scale_fill_manual(values=c("yellow", "#2988E2", "grey50")) +
  scale_y_log10() + 
  coord_flip() + 
  facet_grid(rows = vars(variable), scales="free_y") + 
  theme_classic() + 
  theme(
    text = element_text(size=16)
  )

dev.off()



cat.ub <- as.matrix(read.table("data/processed/abundance_bycategories_ub.txt"))
#merge the RING single + complex data, omit the 'uncharacterized'
cat.ub <- rbind( cat.ub, (cat.ub["RING_Single",] + cat.ub["RING_Complex",]) )
row.names(cat.ub)[9] <- "RING"
cat.ub <- cat.ub[c("E2", "HECT", "RING", "UBOX", "DUB"),]

cat.ub <- as.data.frame(t(cat.ub))
cat.ub$type <- type
cat.ub.m <- melt(cat.ub, id="type")


svg(file = "figures/Figure2/D_total_ub.svg", height = 3.5, width = 4)

ggplot(cat.ub.m) + 
  geom_boxplot(aes(x=type, y=value, fill=type), show.legend = F ) + 
  labs(x="", y="Cumulative expression") + 
  scale_fill_manual(values=c("yellow", "#2988E2", "grey50")) +
  scale_y_log10() + 
  coord_flip() + 
  facet_grid(rows = vars(variable), scales="free_y") + 
  theme_classic() + 
  theme(
     text = element_text(size=16)
  )

dev.off()


  

# Figure 2E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat.all <- as.matrix(read.table("data/processed/cluster_average_norm.txt", header=T))
cat.chap <- as.matrix(read.table("data/processed/abundance_bycategories_chap.txt"))
cat.ub <- as.matrix(read.table("data/processed/abundance_bycategories_ub.txt"))
all( colnames(cat.chap) == colnames(cat.ub) )       # check that cols are corresponding

total.chap <- colSums(cat.chap)
total.ub <- colSums(cat.ub)

total.ratio <- total.chap / total.ub
total.pct <- (total.chap + total.ub) / colSums(cat.all) * 100
df.ratio <- data.frame(total=total.ratio, type=type)
df.pct <- data.frame(total=total.pct, type=type)


plot.ratio <- ggplot(df.ratio, aes(x=total, colour=type)) + 
  stat_ecdf(size=1.5) + 
  labs(x="Chaperones / Ub ligases", y="Fraction") + 
  scale_x_continuous(limits=c(0.45, 0.8)) + 
  scale_color_manual(values=c("yellow", "#2988E2", "grey50")) +
  theme_classic() + 
  theme(
    text = element_text(size=18),
    legend.position = c(0.9, 0.45),
    legend.title = element_blank()
  )


plot.pct <- ggplot(df.pct, aes(x=total, colour=type)) + 
  stat_ecdf(size=1.5, show.legend = F) + 
  labs(x="%PN", y="Fraction") + 
  scale_x_continuous(limits=c(4, 5)) + 
  scale_color_manual(values=c("yellow", "#2988E2", "grey50")) +
  theme_classic() + 
  theme(
    text = element_text(size=18)
)


svg(file = "figures/Figure2/E_pn_ecdf.svg", height = 4.5, width = 3.5)

plot_grid(plot.pct, plot.ratio, labels ="", ncol = 1, align = 'v')

dev.off()




