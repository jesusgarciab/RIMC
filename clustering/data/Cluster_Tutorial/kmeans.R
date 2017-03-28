#kmeans
kclust <- kmeans(scaled_wd,centers = 4,iter.max = 100)

ggplot(comp_dt,aes(Comp.1,Comp.2))+
  geom_point(aes(color = as.factor(kclust$cluster)),size=3)

tunek <- kmeansruns(scaled_wd,krange = 1:10,criterion = "ch")
tunek$bestk #3

tunekw <- kmeansruns(scaled_wd,krange = 1:10,criterion = "asw")
tunekw$bestk #4

