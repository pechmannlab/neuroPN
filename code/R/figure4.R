# Figure 4 neuroPN
# SP | 2022 | sebastian@pechmannlab.net

library(ggplot2)
library(reshape2)
library(ggdendro)
library(cowplot)
library(edgeR)


setwd("M2/neuro/")



# Figure 4B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

coevo <- as.data.frame(read.table("data/processed/coevo/coevoDF.txt", header=T))
coevo <- coevo[c(1,2,4,5,7,8),]   # not plotting chaperone-chaperone and ub-ub interactions, but data is there
coevo$name <- rev(c("Control", "Ribosome", "Ub ligases - synapse", "Ub ligases", "Chaperones - synapse", "Chaperones")) #check that new names match
coevo$name <- factor(coevo$name, levels=c("Control", "Ribosome", "Ub ligases - synapse", "Ub ligases", "Chaperones - synapse", "Chaperones"))
coevo2 <- data.frame(name=coevo[,6],  enrich=(coevo[,2]/coevo[,3]-1)*100 )


svg("figures/Figure4/B_conservation.svg", height=3, width=4)

ggplot(coevo2, aes(x=name, y=enrich)) + 
  geom_col(aes(fill=name), position=position_dodge2(),  show.legend=F ) + 
  scale_fill_manual(values=rev(c("orange", "orange", "purple", "purple", "darkgreen", "grey30"))) + 
  labs(y="Conservation", x="") + 
  scale_y_continuous(breaks=c(0, 25, 50), labels=c(0, 25, 50)) + 
  coord_flip() +
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

dev.off()



# Figure 4C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

degree <- read.table("data/processed/coevo/degree.txt", header=T)
degree$colo <- ifelse(degree$cons > degree$q95, 1 , 0)
degree$colo = factor(degree$colo)


svg("figures/Figure4/C_degree_scatter.svg", height=4, width=3.5)

ggplot(degree) + 
  geom_point(aes(x=rand, y=cons, color=colo ), size=0.8, show.legend = F) + 
  scale_color_manual(values=c("#222222", "#2278AB")) + 
  xlim(c(0, 400)) + 
  ylim(c(0, 500)) + 
  labs(x="Expected degree", y="Conserved degree") + 
  geom_abline(slope=1, intercept=0, size=1.5, color="red") + 
  theme_classic() + 
  theme(
    text = element_text(size=20)
  )

dev.off()




svg("figures/Figure4/C_degree_bar.svg", height=4, width=2)

counts <- data.frame(counts=c(sum(degree$colo==0), sum(degree$colo==1)), cate=c('n.s.', 'sig.')  )
ggplot(counts) + 
  geom_col(aes(x=cate, y=counts, fill=cate), show.legend = T) + 
  scale_fill_manual(values=c("#222222", "#2278AB")) +
  labs(x="", y="Number of genes") + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.title = element_blank(),
    legend.text = element_text(size=14),
    legend.position = c(0.975, 0.9)
  )

dev.off()






# Figure 4D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

degree <- read.table("data/processed/coevo/degree.txt", header=T)
degree$name = row.names(degree)

deg.chap <- degree[degree$cat=="chap",]
deg.chap <- deg.chap[deg.chap$cons > deg.chap$q95, ]
deg.chap$ratio <- deg.chap$cons / deg.chap$rand
tops.chap.cons <- deg.chap[ sort(deg.chap$cons, index.return=T, decreasing=T)$ix, ][1:7,]
tops.chap.ratio <- deg.chap[ sort(deg.chap$ratio, index.return=T, decreasing=T)$ix, ][1:7,]

deg.ub <- degree[degree$cat=="ub",]
deg.ub <- deg.ub[deg.ub$cons > deg.ub$q95, ]
deg.ub$ratio <- deg.ub$cons / deg.ub$rand
tops.ub.cons <- deg.ub[ sort(deg.ub$cons, index.return=T, decreasing=T)$ix, ][1:7,]
tops.ub.ratio <- deg.ub[ sort(deg.ub$ratio, index.return=T, decreasing=T)$ix, ][1:7,]



topsplot <- function(dataIN, AX){
  
  dataIN$name <- factor(dataIN$name, levels=rev(dataIN$name))
  dataIN$conserved <- dataIN$cons - dataIN$rand
  dataIN <- data.frame(name=dataIN$name, conserved=dataIN$conserved, rand=dataIN$rand)
  dataIN <- melt(dataIN, id="name")
  
  ggplot(dataIN) + 
    geom_col(aes(x=name, y=value, fill=variable), show.legend=F) + 
    scale_fill_manual(values=c("#2278AB", "#222222")) +
    scale_y_continuous(breaks=AX, labels=AX) + 
    coord_flip() + 
    labs(x="", y="Degree") + 
    theme_classic() + 
    theme(
      text = element_text(size=16), 
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_blank()
    )
  
}


svg("figures/Figure4/D_tops.svg", height=4, width=5)

plot.chap.cons <- topsplot(tops.chap.cons, c(0, 200, 400))
plot.chap.ratio <- topsplot(tops.chap.ratio, c(0, 100, 200, 300))
plot.ub.cons <- topsplot(tops.ub.cons, c(0, 200, 400))
plot.ub.ratio <- topsplot(tops.ub.ratio, c(0, 100, 200, 300))

plot_grid(plot.chap.cons, plot.chap.ratio, plot.ub.cons, plot.ub.ratio, labels ="", ncol = 2, align = 'h')

dev.off()





# Figure 4E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dens <- data.frame(read.table("data/processed/coevo/density.txt", header=T))
dens$idx <- factor(c("Chaperones", "Ub ligases", "Synapse", "Control"), levels=c("Chaperones", "Ub ligases", "Synapse", "Control") )


rand.chap <- scan("data/processed/coevo/rand_density_chap.txt")
rand.ub <- scan("data/processed/coevo/rand_density_ub.txt")
rand.syn <- scan("data/processed/coevo/rand_density_syn.txt")
rand.rando <- as.matrix(read.table("data/processed/coevo/rand_density_rando.txt"))


df <- rbind(data.frame(rand=rand.chap, cat=rep("Chaperones", length(rand.chap))),
            data.frame(rand=rand.ub, cat=rep("Ub ligases", length(rand.ub))),
            data.frame(rand=rand.syn, cat=rep("Synapse", length(rand.syn))), 
            data.frame(rand=as.numeric(rand.rando), cat=rep("Control", length(as.numeric(rand.rando)))) )


df2 <- data.frame(cat=unique(df$cat), true=dens$true)

svg("figures/Figure4/E_randstar.svg", height=4, width=4)

ggplot(df) + 
  geom_boxplot(aes(x=cat, y=rand, fill=cat), show.legend=F) + 
  geom_point(data=df2, aes(x=cat, y=true), shape=8, size=2, stroke=1.5, color="red") + 
  scale_fill_manual(values=c("orange", "purple", "lightgreen", "grey50")) + 
  labs(x="", y="Network density") + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(angle=40, hjust=1),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )

dev.off()



# Figure 4F ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dens <- data.frame(read.table("data/processed/coevo/density.txt", header=T))
dens$idx <- factor(c("Chaperones", "Ub ligases", "Synapse", "Control"), levels=c("Chaperones", "Ub ligases", "Synapse", "Control") )


svg("figures/Figure4/F_log2ratio.svg", height=4, width=4)

ggplot(dens) + 
  geom_col(aes(x=idx, y=log2(ratio), fill=idx ), show.legend=F) + 
  labs(x="", y="Rel. network density") + 
  scale_fill_manual(values=c("orange", "purple", "grey50", "#222222")) + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(angle=40, hjust=1),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )

dev.off()

