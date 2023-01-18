library(uwot)

human <- as.matrix(read.table("../data/processed/data_merged.txt", header=T))
um <- umap(t(human), metric='correlation', init="random", min_dist = 0.01, n_neighbors=30)

write.table(um, file="../data/processed/umap.txt", quote=F, row.names=F)
