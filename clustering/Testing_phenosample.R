##For tsne
library(RColorBrewer)
library(data.table)
library(fpc)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

setDT(dat[,1:14])

scaled_dat <- scale(dat[,1:14])
    
d <- dist(scaled_dat,method = "euclidean") #distance matrix
h_clust <- hclust(d, method = "ward") #clustering
plot(h_clust)

rect.hclust(h_clust,k=30)

groups <- as.factor(cutree(h_clust,k=30))
dat_and_clusters <- cbind(dat, groups)

ggplot(dat, aes(x = tsne_1, y = tsne_2, color = groups)) + 
  geom_point()

#PCA

pcmp <- princomp(scaled_dat)
pred_pc <- predict(pcmp, newdata=scaled_dat)[,1:2]

comp_dt <- cbind(as.data.table(pred_pc),cluster = as.factor(groups))

ggplot(comp_dt,aes(Comp.1,Comp.2))+
  geom_point(aes(color = cluster),size=3)





##kmeans

kclust <- kmeans(scaled_dat ,centers = 30,iter.max = 100)


ggplot(comp_dt,aes(Comp.1,Comp.2))+
  geom_point(aes(color = as.factor(kclust$cluster)),size=.1)


pamk.best <- pamk(scaled_dat)

tunek <- kmeansruns(scaled_dat,krange = 1:10,criterion = "ch")
tunek$bestk #3

tunekw <- kmeansruns(scaled_dat,krange = 1:10,criterion = "asw")
tunekw$bestk #4

