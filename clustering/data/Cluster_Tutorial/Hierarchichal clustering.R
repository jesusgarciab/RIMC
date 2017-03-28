#set working directory
path <- "data/Cluster_Tutorial"
#setwd(path)


#load data
water_data <- read.table("water-treatment.data.txt",sep = ",",header = F,na.strings = c("?"))
setDT(water_data)
head(water_data)

#lets check missing values
colSums(is.na(water_data))



#impute missing values with median
for(i in colnames(water_data)[!(colnames(water_data) %in% c("V1"))])
set(x = water_data,i = which(is.na(water_data[[i]])), j = i, value = median(water_data[[i]], na.rm = T))

#scale the variables
scaled_wd <- scale(water_data[,-c("V1"),with=F])

#Hierarchical Clustering
d <- dist(scaled_wd,method = "euclidean") #distance matrix
h_clust <- hclust(d, method = "ward") #clustering
plot(h_clust,labels = water_data$V1) #dendrogram


rect.hclust(h_clust,k=4)

#extract clusters
groups <- cutree(h_clust,k=4)
groups

#pca
pcmp <- princomp(scaled_wd)
pred_pc <- predict(pcmp, newdata=scaled_wd)[,1:2]

comp_dt <- cbind(as.data.table(pred_pc),cluster = as.factor(groups), Labels = water_data$V1)

ggplot(comp_dt,aes(Comp.1,Comp.2))+
  geom_point(aes(color = cluster),size=3)


